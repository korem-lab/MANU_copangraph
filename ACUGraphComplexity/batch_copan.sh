#!/bin/bash 
#SBATCH --job-name=cpg
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/log_%x-%j.log
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.005_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.005 -mj 250 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.01_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.01 -mj 750 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.02_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.02 -mj 750 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.03_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.03 -mj 1000 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.04_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.04 -mj 1000 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.05_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.05 -mj 1000 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.06_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.06 -mj 1000 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
/burg/pmg/users/ic2465/copangraph/bin/release/copangraph -g ${REP}_sd_0.1_graph -s /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/${REP}.peext.list -d 0.1 -mj 1000 -ms 100 -o /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/copangraph
