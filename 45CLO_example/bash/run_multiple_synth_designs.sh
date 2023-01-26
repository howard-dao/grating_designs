# run multiple designs

wavelengths=(1270 1580)
#angles=(-20 -15 -10 10 15 20)
#angles=(0)
angles=(-20 20)
etches=('shallow' 'full')
directions=('down')

for wl in ${wavelengths[@]}; do
    for angle in ${angles[@]}; do
        for d in ${directions[@]}; do
            for e in ${etches[@]}; do

                    qsub run_synth_bash_45CLO.sh $wl $angle $e $d 5

            done
        done
    done
done

# # 1270
# qsub run_synth_bash_45CLO.sh 1270 8 shallow down 5
# qsub run_synth_bash_45CLO.sh 1270 -8 shallow down 5
# qsub run_synth_bash_45CLO.sh 1270 8 shallow up 5
# qsub run_synth_bash_45CLO.sh 1270 -8 shallow up 5
# qsub run_synth_bash_45CLO.sh 1270 8 full down 5
# qsub run_synth_bash_45CLO.sh 1270 -8 full down 5
# qsub run_synth_bash_45CLO.sh 1270 8 full up 5
# qsub run_synth_bash_45CLO.sh 1270 -8 full up 5

# qsub run_synth_bash_45RFSOI.sh 1300 10 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 7.5 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 5 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 -5 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 -7.5 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 -10 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 10 5 oxide
# qsub run_synth_bash_45RFSOI.sh 1300 7.5 5 oxide
# qsub run_synth_bash_45RFSOI.sh 1300 5 5 oxide
# qsub run_synth_bash_45RFSOI.sh 1300 -5 5 oxide
# qsub run_synth_bash_45RFSOI.sh 1300 -7.5 5 oxide
# qsub run_synth_bash_45RFSOI.sh 1300 -10 5 oxide
# qsub run_synth_bash_45RFSOI.sh 1300 20 5 air
# qsub run_synth_bash_45RFSOI.sh 1300 -20 5 air

# # 1180
# # qsub run_synth_bash_45RFSOI.sh 1180 15 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 10 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 7.5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 -5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 -7.5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 -10 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 -15 5 air
# # qsub run_synth_bash_45RFSOI.sh 1180 15 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 10 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 7.5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 -5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 -7.5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 -10 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 -15 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1180 20 5 air

# # 1300
# # qsub run_synth_bash_45RFSOI.sh 1300 10 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 7.5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 -5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 -7.5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 -10 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 10 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1300 7.5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1300 5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1300 -5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1300 -7.5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1300 -10 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1300 20 5 air
# # qsub run_synth_bash_45RFSOI.sh 1300 -20 5 air

# qsub run_synth_bash_45RFSOI.sh 1300 -15 5 air up
# qsub run_synth_bash_45RFSOI.sh 1300 15 5 air up
# qsub run_synth_bash_45RFSOI.sh 1300 20 5 air up

# # # 1550
# # qsub run_synth_bash_45RFSOI.sh 1550 10 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 7.5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 -5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 -7.5 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 -10 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 10 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1550 7.5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1550 5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1550 -5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1550 -7.5 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1550 -10 5 oxide
# # qsub run_synth_bash_45RFSOI.sh 1550 20 5 air
# # qsub run_synth_bash_45RFSOI.sh 1550 -20 5 air

# qsub run_synth_bash_45RFSOI.sh 1550 -15 5 air up
# qsub run_synth_bash_45RFSOI.sh 1550 15 5 air up
# qsub run_synth_bash_45RFSOI.sh 1550 20 5 air up