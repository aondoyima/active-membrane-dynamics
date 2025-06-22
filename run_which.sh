# Usage: bash run_which.sh -script select_script.py

# Load parameter file
source params.txt

# Parse command line argument
while [[ $# -gt 0 ]]; do
    case "$1" in
        -script)
            script="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$script" ]]; then
    echo "Error: No script specified. Use -script <script_name.py>"
    exit 1
fi

# Loop over parameter combinations
for L in $L_list; do
    for zeta in $zeta_list; do
        for kappa in $kappa_list; do
            for sigma in $sigma_list; do
                for K in $K_list; do
                    for chi in $chi_list; do
                        for alpha in $alpha_list; do
                            for s1 in $s1_list; do
                                s2=$s1  # hardcoded

                                if [[ "$script" == "pbc_main.py" ]]; then
                                    python pbc_main.py \
                                        -a $a -b $b -cp $chiprime -K $K \
                                        -s $sigma -k $kappa -chi $chi -z $zeta \
                                        -alpha $alpha -s0 $s0 -s1 $s1 -s2 $s2 \
                                        -psi1bar $psi1bar -psi2bar $psi2bar \
                                        -syslength $L -initnoise $initnoise \
                                        -vc $vol_conserve -ds $data_save \
                                        -ad $adapt_step -vn $vary_noise \
                                        -res $resolution

                                elif [[ "$script" == "pbc_mov.py" ]]; then
                                    python pbc_mov.py \
                                        -a $a -b $b -cp $chiprime -K $K \
                                        -s $sigma -k $kappa -chi $chi -z $zeta \
                                        -alpha $alpha -s0 $s0 -s1 $s1 -s2 $s2 \
                                        -psi1bar $psi1bar -psi2bar $psi2bar \
                                        -syslength $L -initnoise $initnoise \
                                        -vc $vol_conserve -ds $data_save \
                                        -ad $adapt_step -vn $vary_noise \
                                        -mov $mov -kym $kym -amp $kym_amp \
                                        -struct $struct -t1 $t1 -t2 $t2 -tread $t_read

                                elif [[ "$script" == "analysis_sym_reg.py" ]]; then
                                    python analysis_sym_reg.py \
                                        -a $a -b $b -cp $chiprime -K $K \
                                        -s $sigma -k $kappa -chi $chi -z $zeta \
                                        -alpha $alpha -s0 $s0 -s1 $s1 -s2 $s2 \
                                        -psi1bar $psi1bar -psi2bar $psi2bar \
                                        -syslength $L -initnoise $initnoise \
                                        -vc $vol_conserve -ds $data_save \
                                        -ad $adapt_step -vn $vary_noise

                                else
                                    echo "Error: Unknown script '$script'"
                                    exit 1
                                fi

                            done
                        done
                    done
                done
            done
        done
    done
done