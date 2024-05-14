import sys
import random
import evaluation.utils.parse_seq as ps

random.seed(101)
DNA = 'ATCG'


def gen_random_seq(l):
    return ''.join(random.choice(DNA) for _ in range(l))


def num_same(s1, s2):
    return sum(a == b for a, b in zip(s1, s2))

def mutate(seq):
    per_mut = random.randint(75, 100)
    mut = ''.join(
        c if random.randint(1, 100) <= per_mut else random.choice(list(set(DNA) - set(c))) for c in seq
    )

    return num_same(seq, mut) / len(mut), mut


def generate_random_batch(batch_name, l=10000, n=1000):

    seed = ps.Fasta('seed', gen_random_seq(l))
    batch = [seed]
    for i in range(n):
        perc_id, mut = mutate(seed.seq)
        batch.append(ps.Fasta(f'{batch_name}_mutant_{i}_{perc_id}', mut))
    with open(f'{batch_name}.fasta', 'w') as f:
        for fa in batch:
            fa.write(f)


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('Usage: <exe> <batch_base> <num_batches>')
        sys.exit(-1)

    num_batches = int(sys.argv[2])
    batch_base = sys.argv[1]

    for i in range(num_batches):
        print('constructing batch', i)
        generate_random_batch(f'batch_{i}')



