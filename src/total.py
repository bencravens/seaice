import grab
import plot

myvars = ['sithick','sispeed','sihc','siflswdbot','siflsendupbot','siflcondtop','siflcondbot','fsurf_ai','fhocn_ai','ardg','dardg1dt','opening','snoice']
control = 'u-at053'
models = ['u-au866','u-au872','u-au874','u-av231']
months=[2,9]
for var in myvars:
    for month in months:
        plot.month_map_mean_main(control,month,var,"/home/ben/Documents/summer2019/plotlims",False)
        for model in models:
            plot.month_map_anom_main(model,month,var,"/home/ben/Documents/summer2019/plotlims",False)
