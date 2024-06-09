#! /bin/bash

# Description: This script solves the model and computes shock elasticities for the first specification
# in Section 4.5, where ambiguity is limited to the two slope parameters.

#### Parameter Configuration Arrays: each array contains the values of the parameters for the different scenarios
Deltaarray=(1.0)
rhoarray=(1.0)
alphaarray=(0.0922)
qarray=(0.2)
gammaarray=(1.0 3.0 1.0 4.0)
twoparameterarray=(1 1 0 0)

#### Script Names: "solve_name" solves the model; "elasticity_name" computes shock elasticities;
solve_name="main_onecapital.jl"
elasticity_name="main_pde_shock_elasticity.py"

for Delta in ${Deltaarray[@]}; do
    for index in ${!rhoarray[@]}; do
        rho=${rhoarray[$index]}
        alpha=${alphaarray[$index]}
        for index2 in ${!gammaarray[@]}; do
            gamma=${gammaarray[$index2]}
            twoparameter=${twoparameterarray[$index2]}
                for q in "${qarray[@]}"; do

                    action_name="onecap_model_with_structure_ambiguity"

                    mkdir -p ./job-outs/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/

                    if [ -f ./bash/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.sh ]; then
                        rm ./bash/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.sh
                    fi

                    mkdir -p ./bash/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/

                    touch ./bash/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.sh

                    tee -a ./bash/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.sh <<EOF
#!/bin/bash

#SBATCH --account=ssd
#SBATCH --job-name=Delta_${Delta}_twoparameter_${twoparameter}_gamma_${gamma}_q_${q}_rho_${rho}_${action_name}
#SBATCH --output=./job-outs/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.out
#SBATCH --error=./job-outs/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.err
#SBATCH --time=0-23:00:00
#SBATCH --partition=ssd
#SBATCH --qos=ssd
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

module load julia/1.7.3
module load python/anaconda-2020.11

echo "\$SLURM_JOB_NAME"

echo "Program starts \$(date)"
start_time=\$(date +%s)

srun julia ./src/$solve_name  --Delta ${Delta} --alpha ${alpha} --gamma ${gamma} --rho ${rho} --q ${q} --action_name ${action_name}  --twoparameter ${twoparameter} 
python3 ./src/$elasticity_name  --Delta ${Delta} --alpha ${alpha} --gamma ${gamma} --rho ${rho} --q ${q} --action_name ${action_name} --twoparameter ${twoparameter} 

echo "Program ends \$(date)"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))

eval "echo Elapsed time: \$(date -ud "@\$elapsed" +'\$((%s/3600/24)) days %H hr %M min %S sec')"

EOF
                count=$(($count + 1))
                sbatch ./bash/${action_name}/Delta_${Delta}_twoparameter_${twoparameter}/rho_${rho}_gamma_${gamma}_q_${q}_alpha_${alpha}/run.sh
            done
        done
    done
done