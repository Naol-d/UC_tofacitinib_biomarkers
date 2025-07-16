#!/bin/bash -l
#SBATCH --job-name=firefox-gui
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out

readonly PORT=5901
readonly SIF_PATH="/scratch/prj/bmb_tofacitinib/research_project/singularity_images/centosvnc.sif"
readonly WORKSPACE_PATH="/scratch/prj/bmb_tofacitinib/research_project/deconvolution/cibersortx/scripts/browser/internet_access/log/nginx"
readonly IPADDRESS=$(hostname -I | tr ' ' '\n' | grep '10.211.4.')

echo "IP address: $IPADDRESS"


mkdir -p "/scratch/prj/bmb_tofacitinib/research_project/deconvolution_cibersortx/scripts/browser/internet_access/log/nginx"

echo ""
echo "SSH Tunnel..."
echo "ssh -m hmac-sha2-512 -NL ${PORT}:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk"
echo ""
echo "Login details"
echo "Password = Passw0rd"


singularity exec \
  --bind ${WORKSPACE_PATH}:${WORKSPACE_PATH} \
  $SIF_PATH \
  dbus-launch /opt/vnc_startup.sh VNC_PORT:5901 VNC_PW:Passw0rd
