#/bin/bash
#for (( i=36; i<100; i++ )); do DLkl_daqrcnp /data4/cdcSpring8/root run_0000${i}_built.root 1000; done
#for (( i=100; i<126; i++ )); do DLkl_daqrcnp /data4/cdcSpring8/root run_000${i}_built.root 1000; done
#for (( i=1000; i<1079; i++ )); do DLkl_daqrcnp /data4/cdcSpring8/root run_00${i}_built.root 1000; done

for (( i=2001; i<2071; i++))
do
#	DLkl_daqrcnp /data4/cdcSpring8/root-preview run_00${i}_built.root 1000
	scp hiroki@192.168.2.30:/data4/cdcSpring8/root/run_00${i}_built.root root/
done
