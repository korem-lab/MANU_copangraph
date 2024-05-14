import subprocess
import pandas as pd

class NodeBalancer:
    def __init__(self, limit=96):
        self.nodes = {n:0 for n in list(get_node_info().index)}
        self.limit = limit
        self.user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()

        self.reset()

    def get_node(self, num_cpus, exclude=None, include=None):
        avail_nodes = find_nodes(num_cpus).sample(frac=1)
        for n in avail_nodes.index:
            if exclude is not None and n in exclude:
                print('Skipping %s because it is excluded...'%n)
                continue
            if include is not None and n not in include:
                print('Skipping %s because it is not included...'%n)
                continue
            if (self.nodes[n] <= self.limit - num_cpus):
                self.nodes[n] += num_cpus
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
                    new_node = self.get_node(requested_cpus,exclude=[requested_node,'m012','m013'])
                    print('Updated node for %s (%s) from %s to %s'%(job_name,job_id,requested_node,new_node))
                    self.nodes[requested_node] -= requested_cpus

                    subprocess.run(["scontrol","update",f"JobID={job_id}",f"NODELIST={new_node}"]) 
            except:
                continue

def get_node_info():
    pestat_out = subprocess.run(['sinfo', '-o', '%n %e %C'], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    df = []

    for line in pestat_out[1:]:
        line_split = line.split(' ')
        try:
            node = line_split[0]
            mem = round(int(line_split[1])/1000,2)
            used = int(line_split[2].split('/')[0])
            available = int(line_split[2].split('/')[1])
            df.append({'node':node,'memory (GB)':mem,'total':used+available,'used':used,'available':available})
        except:
            continue
        
    df = pd.DataFrame(df).set_index('node')
    return df

def find_nodes(num_cpus):
    node_info = get_node_info()
    node_info = node_info.loc[node_info['available']>=num_cpus]
    return node_info

def get_jobs():
    user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    jobs = subprocess.run(['squeue','-o','%.30j','-u', user], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    return [l.strip() for l in jobs[1:]]

def get_job_info():
    user = subprocess.run(['whoami'], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    jobs = subprocess.run(['squeue','-o','%.18i %.30j','-u', user], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
    return [l.strip() for l in jobs[1:]]

if __name__ == "__main__":
    nodeinfo = get_node_info().sort_values('available',ascending=False)
    df = pd.concat([nodeinfo,NodeBalancer().get_usage()],axis=1)
    df = df.loc[sorted(df.index)]
    df.loc['-----']=pd.Series()
    df.loc['TOTAL']=df.sum()

    print(df.fillna('-----'))
