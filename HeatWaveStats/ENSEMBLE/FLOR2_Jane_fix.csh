dmget /archive/jwb/STEP/ENS2_FIX/RAW/atmos_daily_1941-2050_t_ref*ens$1.nc
cp /archive/jwb/STEP/ENS2_FIX/RAW/atmos_daily_1941-2050_t_ref*ens$1.nc $TMPDIR
source /home/Jane.Baldwin/anaconda2/bin/activate /home/Jane.Baldwin/anaconda2/envs/mypy2

python /home/Jane.Baldwin/PYTHON/STEP/ehfheatwaves_compound_3def_fix7-18-17.py -x atmos_daily_1941-2050_t_ref_max_ens$1.nc --vnamex=t_ref_max -n atmos_daily_1941-2050_t_ref_min_ens$1.nc --vnamen=t_ref_min --vnamet=time -m /ptmp/jwb/atmos_daily.static.nc --vnamem=land_mask --t90pc --base=1981-2010 -d FLOR -e $1 -p 90
cp tx90*r$1* tn90*r$1* EHF90*r$1* /archive/jwb/STEP/ENS2_FIX/PROC/

rm -f /ptmp/jwb/atmos_daily_1941-2050_t_ref*ens$1.nc















