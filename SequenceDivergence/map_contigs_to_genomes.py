import os
import sys
import subprocess
import pickle
import re
import yaml
import pandas as pd
import numpy as np
import parse_seq as ps


class SubProcessTools:
    makeblastdb = 'makeblastdb'
    blastn = 'blastn'


class GroundTruth:
    """Constructs a ground truth matrix (genomes by contigs) from an input blast file. """

    
    def __init__(self, yaml_config_path):


        # Parse config file
        with open(yaml_config_path) as f:
            self.config = yaml.load(f, Loader=yaml.SafeLoader)

        self.min_contig_len = self.config['min_contig_len']
        self.min_pident = self.config['min_pident']
        self.min_len_match = self.config['min_len_match']
        self.contigs = self.config['contigs']
        self.reference_genomes = self.config['reference_genomes']
        self.dataset_name = self.config['dataset_name']
        self.out_dir = self.config['out_dir']
        self.blast_dir = self.config['blast_dir']
        self.asm_tool = self.config['asm_tool']

        # Get the blast file. If it doesn't exist, then run blast (i.e make it)
        self.blast = self.get_blast_file()

        self.genome_names = pd.Series(sorted(self.get_genome_names()))

        # Construct blast hit matrix.
        self.filter_blast_hits()

        # Get names of remaining contigs.
        # Contigs are in fastg format. The header has structure >EDGE_100_blah_blah_blah
        # We capture the integer after EDGE which then becomes the contig name / id
        # TODO: Rather than filtering out contigs, Add a binary, "modelled" column, which is
        # TODO: True if the contig is to be modelled, and false otherwise
        self.contig_info = self.get_contig_info_from_blast(self.asm_tool)

        # Construct the ground truth.
        self.ground_truth = self.construct_ground_truth()

        #self.sp_lvl_gt = self.construct_species_level_ground_truth()

        # Write out data
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        self.blast.to_pickle(os.path.join(self.out_dir, f'{self.dataset_name}_blast.pkl'))
        self.ground_truth.to_pickle(os.path.join(self.out_dir, f'{self.dataset_name}_ground_truth.pkl'))
        #self.sp_lvl_gt.to_pickle(os.path.join(self.out_dir, 'sp_level_gt.pkl'))
        self.contig_info.to_pickle(os.path.join(self.out_dir, f'{self.dataset_name}_contig_info.pkl'))
        genome_file = os.path.join(self.out_dir, f'{self.dataset_name}_genome_map.pkl')
        self.genome_names = pd.Series(self.genome_names)
        self.genome_names.to_pickle(genome_file)

    def construct_species_level_ground_truth(self):
        self.taxonomy_map = pd.read_pickle(self.taxonomy_map)
        self.taxonomy_map = self.taxonomy_map.drop_duplicates()
        self.taxonomy_map.species = self.taxonomy_map.species.astype(int)
        self.taxonomy_map = self.taxonomy_map.sort_values(by='filename')
        self.taxonomy_map.loc[:, 'genomes'] = self.ground_truth.index
        species = self.taxonomy_map.species.unique()
        species.sort()

        # Species table is going to be Species X contigs
        sp_gt = pd.DataFrame(
            np.full((len(species), len(self.ground_truth.columns)), False),
            index=species,
            columns=self.ground_truth.columns
        )

        # iterate through each species. Grab all genomes assigned
        # to species sp and take the logical OR of their ground truths.
        # This computes a boolean vector of length contigs, describing
        # which contigs are present in species sp.
        for sp in species:

            # Get genomes associated with sp
            genomes_in_sp = self.taxonomy_map.loc[
                self.taxonomy_map.species == sp, :
            ].genomes

            # OR all genome ground truth vectors assigned to sp
            sp_gt.loc[sp, :] = self.ground_truth.loc[genomes_in_sp, :].any()

        return sp_gt

    def get_contig_info_from_blast(self, contig_source):

        # Extract a unique integer id for every contig.
        # Since the format of the fasta differs by algorithm
        # a different regex will be needed per algorithm
        if contig_source == 'ms':
            contig_ids = [
                int(re.findall('NODE_([0-9]+)_', s)[0]) for s in self.blast.contig
            ]
        elif contig_source == 'mh':
            contig_ids = [
                int(re.findall('k[0-9]+_([0-9]+)$', s)[0]) for s in self.blast.contig
            ]
        elif contig_source == 'ms_graph':
            contig_ids = [
                int(re.findall('EDGE_([0-9]+)_', s)[0]) for s in self.blast.contig
            ]
        # Constructing a dataframe from a list of lists will set each list
        # to a row. I want each list to be a column, so transpose the dataframe
        info = pd.DataFrame(
            [
                self.blast.contig,
                self.blast.contig_length
            ]
        ).T

        # Set the correct column value, and then drop duplicate rows.
        # Duplicates are present because of the blast file - the same contig
        # can map to multiple locations, so the same contig can have multiple entries
        info.columns = ['raw_contig_id', 'contig_length']
        info.index = contig_ids
        info.index.name = 'contig_id'
        return info.drop_duplicates()

    def get_genome_names(self):
        return [
            seq.hdr.split()[0] for seq in ps.parse(open(self.reference_genomes), ps.Fasta)
        ]

    def get_blast_file(self, contig_source='metaspades'):
        """Blasts contigs against reference genomes."""
        blast_results = os.path.join(self.config['blast_dir'], f'{self.dataset_name}_blast_results.txt')
        if not os.path.exists(blast_results):
            # make the blast database
            blast_dir = self.blast_dir
            blast_db = os.path.join(blast_dir, self.dataset_name)
            if not os.path.exists(blast_dir):
                os.mkdir(blast_dir)
            cmd = f'{SubProcessTools.makeblastdb} -dbtype nucl ' + \
                  f'-in {self.reference_genomes} -out {blast_db}'
            subprocess.call(cmd, shell=True)

            # run blast
            cmd = f'{SubProcessTools.blastn} -db {blast_db} ' + \
                  f'-outfmt "6 qacc sacc evalue qstart qend qlen sstart send slen pident nident" ' + \
                  f'-query {self.contigs} -out {blast_results}'
            subprocess.call(cmd, shell=True)

        # Load blast result into pandas dataframe
        blast = pd.read_csv(
            blast_results,
            delimiter='\t',
            header=None
        )
        blast.columns = [
            'contig', 'genome', 'evalue', 'contig_start', 'contig_end',
            'contig_length', 'genome_start', 'genome_end', 'genome_length',
            'pident', 'nident'
        ]

        return blast

    def filter_blast_hits(self):
        """
        Filters hits to those assumed to be the true location
        for which a contig derived (contigs can map to multiple genomes).
        """

        filter = self.blast.contig_length >= self.min_contig_len
        # The contig must be a near identical match to the genome
        filter &= self.blast.pident >= self.min_pident
        # The length of the match must be nearly the entire length of the contig
        # TODO: Come up with a function that match must be identical within n characters.
        # TODO: This avoids the issue of very long contigs being allowed a large amount of unmatching sequence.
        # TODO: n needs to decrease for short contigs.
        filter &= self.blast.nident / self.blast.contig_length * 100 >= self.min_len_match
        self.blast = self.blast.loc[filter, :]

    def construct_ground_truth(self):
        """
        Converts the blast hits into a contigs X genomes matrix,
        where entry i,j is true if contig i mapped to genome j.
        """

        # Initialize ground truth matrix
        G, C = len(self.genome_names), len(self.contig_info)
        gt = pd.DataFrame(
            np.full((G, C), False),
            index=self.genome_names,
            columns=self.contig_info.raw_contig_id
        )

        # Loop through each blast hit and find the genome to which the contig
        # mapped, and set the entry in the ground truth matrix to true.
        for idx in self.blast.index:
            gt.loc[self.blast.loc[idx, 'genome'], self.blast.loc[idx, 'contig']] = True

        # Ensure each contig is mapped to at least one genome before returning.
        assert gt.any(axis=0).all(), \
            Exception(f'The following contigs were not mapped to any genome {gt.index[~gt.any()]}')

        # Replace genome name with integer value (the map between the two is stored)
        gt.index = self.genome_names.index
        gt.columns = self.contig_info.index
        return gt


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: <exe> config')
        sys.exit()

    _, config = sys.argv
    print('Blasting...', flush=True)
    gt = GroundTruth(config)
    tt = pd.read_csv(gt.config["taxonomic_table"], index_col=0)
    f2c = gt.ground_truth
    f2c.columns = gt.contig_info.raw_contig_id
    f2c.index = gt.genome_names
    num_genomes = tt.shape[0]
    num_contigs = f2c.shape[1]

    print('Translating fragments to genomes...', flush=True)
    g2c = pd.DataFrame(np.full((num_genomes, num_contigs), False), index=tt.file_name, columns=f2c.columns)
    g2c.index.name = 'genome'
    frags_list = list()
    for fl in tt.file_name:
        frags = {e.hdr.split()[0] for e in ps.parse(open(os.path.join(gt.config["genome_dir"], fl)), ps.Fasta)}
        frags_list.append((fl, frags))
        for frag in frags:
            g2c.loc[fl, :] |= f2c.loc[frag, :]
    g2c.columns = gt.contig_info.raw_contig_id
    g2c.to_csv(os.path.join(gt.config['out_dir'], f'{gt.config["dataset_name"]}_c2g.csv'))
    pickle.dump(frags_list, open(os.path.join(gt.config['out_dir'], f'{gt.config["dataset_name"]}_frag_list.pkl'), 'wb'))




