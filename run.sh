export HOME1=/opt/ofsdev/users/azhang/SKILL_OFS
cd $HOME1
cp my_parameters.ctl.cbofs.wl my_parameters.ctl
/opt/ofsdev/users/azhang/SKILL_OFS/scripts/SKILLSTEPS.sh
echo "running SA for temperature"
cd $HOME1
cp my_parameters.ctl.cbofs.temp my_parameters.ctl
/opt/ofsdev/users/azhang/SKILL_OFS/scripts/SKILLSTEPS.sh
echo "running SA for salinity"
cd $HOME1
cp my_parameters.ctl.cbofs.salt my_parameters.ctl
/opt/ofsdev/users/azhang/SKILL_OFS/scripts/SKILLSTEPS.sh

