import Utils.AssemblyParser as AP
import sys

if __name__ == '__main__':
    graph = sys.argv[1]
    asm = AP.Assembly('graph', graph)
    asm.get_max_node_lens().to_csv(graph.replace('.gfa', '').replace('.fastg', '') + '.max_len.csv')
