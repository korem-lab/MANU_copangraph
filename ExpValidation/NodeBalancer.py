import subprocess
import pandas as pd
import subprocess

from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

class NodeBalancer:
    def __init__(self, limit=96, include=None, exclude=None, force=False):
        self.nodes = {n:0 for n in list(get_node_info().index)}
        self.limit = limit
        self.user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
        self.include = include
        self.exclude = exclude
        self.force = force
        self.reset()

    def get_node(self, num_cpus, exclude=None, include=None, force=None, storage_space=None):
        exclude = exclude if exclude is not None else self.exclude
        include = include if include is not None else self.include
        force = force if force is not None else self.force

        if force and include:
            return include[0]
        
        avail_nodes = find_nodes(-1 if force else num_cpus).sort_values('available',ascending=False)
        print('%s available nodes!'%len(avail_nodes))

        if storage_space != None:
            print(f"Checking nodes for storage usage less than {storage_space}%")
            if not 'diskinfo' in self.__dict__:
                self.diskinfo = get_disk_info()
            nodes = self.diskinfo[self.diskinfo['usage'].str.split('%').str[0].astype(float)<storage_space].index
            avail_nodes = avail_nodes.loc[avail_nodes.index.intersection(nodes)]

        for n in avail_nodes.sample(frac=1).index:
            if exclude is not None and n in exclude:
                print('Skipping %s because it is excluded...'%n)
                continue
            if include is not None and n not in include:
                print('Skipping %s because it is not in %s'%(n,include))
                continue
            if force or (self.nodes[n] <= self.limit - num_cpus):
                self.nodes[n] += num_cpus
                print('Returning node %s'%n)
                return n
        raise Exception('No available nodes left! Check your usage...')

    def reset(self):
        self.nodes = {n:0 for n in list(get_node_info().index)}

        pestat_out = subprocess.run(['squeue','-o','%5N %.5C','-u', self.user], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')

        for line in pestat_out[1:]:
            line_split = line.split(' ')
            try:
                node = str(line_split[0])
                cpus = int(line_split[-1])

                self.nodes[node] += cpus
            except:
                continue
    
    def get_usage(self):
        return pd.DataFrame({self.user:self.nodes})
    
    def update_stale_jobs(self):
        user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
        jobs = subprocess.run(['squeue','-o','%.18i %.30j %.4n %.3C %.2t','-u', user], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
        for j in jobs[1:]:
            job_split = [x for x in j.split(' ') if x != '']
            try:
                job_id = job_split[0]
                job_name = job_split[1]
                requested_node = job_split[2]
                requested_cpus = int(job_split[3])
                status = job_split[4]

                if status == 'PD':
                    new_node = self.get_node(requested_cpus,exclude=[requested_node,'m012','m013','m014','m015'])
                    print('Updated node for %s (%s) from %s to %s'%(job_name,job_id,requested_node,new_node))
                    self.nodes[requested_node] -= requested_cpus
                    subprocess.run(["scontrol","update",f"JobID={job_id}",f"NODELIST={new_node}"]) 
            except:
                continue

def get_disk_info():
    nodes = get_node_info()
    storage_info = []
    for node in nodes[nodes['state']=='mix'].index:
        print('Checking node %s'%node)
        with suppress_stdout():
            try:
                df_out = subprocess.run(['srun','--account=pmg','--time=0:0:30','--nodelist='+node,'df','-h'], stdout=subprocess.PIPE).stdout.decode('utf-8')
                line = [x for x in [l for l in df_out.split('\n') if '/pmglocal' in l][0].split(' ') if x != '']
                percent = line[-2]
                avail = line[-3]
                used = line[-4]
                storage_info.append({'node':node,'used':used,'avail':avail,'usage':percent})
            except:
                continue
    df = pd.DataFrame(storage_info).set_index('node').sort_values(by='usage')
    df.index.name = ''
    return df

def get_node_info():
    pestat_out = subprocess.run(['sinfo', '-o', '%n %e %C %t'], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    df = []

    for line in pestat_out[1:]:
        line_split = line.split(' ')
        try:
            node = line_split[0]
            mem = round(int(line_split[1])/1000,2)
            used = int(line_split[2].split('/')[0])
            available = int(line_split[2].split('/')[1])
            state = line_split[3]
            df.append({'node':node,'memory (GB)':mem,'total':used+available,'used':used,'available':available,'state':state})
        except:
            continue
        
    df = pd.DataFrame(df).set_index('node')
    return df

def find_nodes(num_cpus):
    node_info = get_node_info()
    node_info = node_info.loc[(node_info['available']>=num_cpus)&(node_info['state']=='mix')]
    return node_info

def get_jobs():
    user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    jobs = subprocess.run(['squeue','-o','%.99j','-u', user], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    return [l.strip() for l in jobs[1:]]

def get_job_info():
    user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    jobs = subprocess.run(['squeue','-o','%.18i %.30j','-u', user], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    return [l.strip() for l in jobs[1:]]

if __name__ == "__main__":
    nodeinfo = get_node_info().sort_values('available',ascending=False)
    df = pd.concat([nodeinfo.drop('state',axis=1),NodeBalancer().get_usage()],axis=1)
    df = df.loc[sorted(df.index)]
    df.loc['-----']=pd.Series()
    df.loc['TOTAL']=df.sum()
    df = pd.concat([df,nodeinfo['state']],axis=1)

    print(df.fillna('-----'))
