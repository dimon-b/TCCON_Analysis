# -*- coding: utf-8 -*-
"""
Created on 18/12/21 08:10
Project    TCCON analysis
@author:   Dmitry
    Work with TCCON observation data
"""

import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import pickle
# import pytz
from os import listdir


class do_all():

    #
    # =================================== init
    def __init__(self):
        self.wrk_dir = "D:/OneDrive - ees.hokudai.ac.jp/HokkaidoWork/TCOON/"
        self.inp_dir = self.wrk_dir + '2_tc-dset/1_tc-orig/'
        self.mdm_dir = self.wrk_dir + '2_tc-dset/2_ts-pickle/'
        self.res_dir = self.wrk_dir + '3_tc-plot/'

    #
    # =================================== run
    def run(self):
        print('\n\tWorking dir : ', self.wrk_dir)

        # === pd expand option
        pd.set_option('expand_frame_repr', False)

        # === read original files
        self.get_row()

    # =================================== run Prophet
    def twr_pht(self):

        print('\n <<<<< Obs df analysis with prothet >>>>>')

        # === ph changepoint_prior_scale
        cps = 0.01  # 0.005 # less to smoth
        cps = 0.01  # 0.02  more variation for IGR

        # === yearly_seasonality
        yrs = 5

        # === trend
        tr = 0
        if tr:
            self.pht_trend(cps)

        # === seasonal
        se = 0
        if se:
            self.pht_seas(cps, yrs)

        # === weekly
        wk = 1
        if wk:
            self.pht_weekly(cps)
            # self.pht_weekly_s(cps)

        # === hourly amplitude
        hr = 0
        if hr:
            amp_f = open(self.pdir + 'ph/amp.txt', "w")
            self.pht_hourly_s4(cps, amp_f)
            amp_f.close()

        # === hourly
        hr = 1
        if hr:
            self.pht_hourly(cps)
            # self.pht_hourly_sm12(cps)

        # --- all
        al = 0
        if al:
            self.prophet(cps)

    # =================================== run Prophet trend
    def pht_trend(self, cps):

        # --- resample
        rsm_time = 'w'

        # --- loop sites and tracers
        for itr in range(0, len(self.tracer), 1):

            # --- fig
            fig = plt.figure(figsize=(9, 6))
            plt.rc('font', family='serif')
            font_size = 16
            mpl.rcParams.update({'font.size': font_size})

            for ist in range(0, len(self.sites), 1):
                print('\n\n\t' + self.sites[ist][0])

                # --- cps
                cps_n = cps
                if self.sites[ist][0] == 'DEM':
                    cps_n = 0.004

                # --- read df
                df = self.get_1_df(self.sites[ist][0], self.tracer[itr])

                # --- preproc
                df, cols = self.preproc(df, self.tracer[itr])

                if rsm_time != 'n':
                    df = df.resample(rsm_time).mean()

                df['ds'] = df.index
                if self.tracer[itr] == 'CH4':
                    df.rename(columns={'CH4_L': 'y'}, inplace=True)
                    ylab = 'CH$_4$, ppb'
                    ylims = [1900, 2100]
                else:
                    df.rename(columns={'CO2_L': 'y'}, inplace=True)
                    ylab = 'CO$_2$, ppm'
                    ylims = [380, 410]

                # --- ph model
                ph = Prophet(changepoint_prior_scale=cps_n)
                ph.fit(df)
                future = ph.make_future_dataframe(periods=0)
                future.tail()
                fc = ph.predict(future)

                # --- plot
                plt.plot(fc['ds'].dt.to_pydatetime(), fc['trend'], linewidth=2, label=self.sites[ist][0])

            axes = plt.gca()
            axes.set_xlim(pd.Timestamp('2004-11-01'), pd.Timestamp('2015-03-01'))
            axes.set_ylim(ylims)
            axes.set_xlabel('Year')
            axes.set_ylabel(ylab)

            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
            plt.legend(loc="upper left", ncol=2)
            plt.grid(True)

            acc_saveplot.save_plot(self.pdir + 'trend_' + self.tracer[itr], ext="png", close=True, verbose=False)

    # --- prophet seasonal



    #
    # =================================== get picle data
    def get_picles(self):

        # === atrib
        with open(self.mdm_dir + 'all_atrib', 'rb') as fp:
            all_atrib = pickle.load(fp)

        # === dataframe
        all_ds = pd.read_pickle(self.mdm_dir + 'all_ds')

    #
    # =================================== get data from nc files
    def get_row(self):
        #
        # =============================== dir ls
        def dir_ls(inp_dir):
            inp_files = listdir(inp_dir)
            n_files = len(inp_files)
            print('\tNumber of input files == ', n_files)
            return inp_files, n_files

        # =============================== datetime
        def make_date(year, day, hour):  # , tzdata=pytz.UTC):
            dt = datetime.datetime(year, 1, 1)  # , tzdata=tzdata)
            return dt + datetime.timedelta(days=int(day), hours=float(hour))

        # === get list of files
        self.inp_files, self.n_files = dir_ls(self.inp_dir)

        all_atrib = []

        # === loop sites
        for i in range(0, self.n_files, 1):
            n_file = self.inp_dir + self.inp_files[i]
            nc = netCDF4.Dataset(n_file)

            # === print file info
            info = 0
            if info:
                print('\n\tAttributes: ')
                print(nc.ncattrs())
                print('\n\tKeys: ')
                print(nc.variables.keys())

            # === get attributes
            atrib = [nc.getncattr('longName'), nc.getncattr('Location'),
                     round(float(nc.variables['long_deg'][0]), 2),
                     round(float(nc.variables['lat_deg'][0]), 2)]
            all_atrib.append(atrib)

            # === get variables
            year = nc.variables['year'][:]
            day = nc.variables['day'][:]
            hour = nc.variables['hour'][:]
            xco2_ppm = nc.variables['xco2_ppm'][:]

            # === datatime from year, day, hour
            dts = [make_date(yr, dy, hr) for yr, dy, hr in zip(year, day, hour)]

            # === print attributes info
            if info:
                print('\n\tSite attributes: ')
                print(atrib)

            # === Create DataFrame
            ds = pd.DataFrame({nc.getncattr('longName'): xco2_ppm}, index=dts)
            ds = ds.resample(rule='1H').mean()

            # === merge
            if i == 0:
                res_ds = ds
            else:
                # === ds header
                if info:
                    print(ds.head())
                    print(res_ds.head())

                # === concat
                res_ds = pd.concat([res_ds, ds], axis=1)

        # === result head
        print(res_ds.head())

        # === all atrib
        print(all_atrib)

        # === test plot
        res_ds.plot()
        self.save_plot(self.res_dir + 'test_ds.png')

        # === save to file
        res_ds.to_pickle(self.mdm_dir + 'all_ds')
        with open(self.mdm_dir + 'all_atrib', 'wb') as fp:
            pickle.dump(all_atrib, fp)

    #
    # ======================================= save plot
    def save_plot(self, path, ext='png', close=True, verbose=False):
        import matplotlib.pyplot as plt
        import os

        # --- Extract the directory and filename from the given path
        directory = os.path.split(path)[0]
        filename = "%s.%s"%(os.path.split(path)[1], ext)
        if directory == '':
            directory = '.'

        # --- If the directory does not exist, create it
        if not os.path.exists(directory):
            os.makedirs(directory)

        # --- The final path to save to
        savepath = os.path.join(directory, filename)

        if verbose:
            print('\n\tSaving figure to           : ', savepath)

        # --- Actually save the figure
        plt.savefig(savepath, dpi=300, bbox_inches='tight')

        # --- options
        if close:
            plt.close()
        if verbose:
            print("plot_save - OK")


#
#
# ======================================= main
if __name__ == '__main__':

    #
    # =================================== estimate time
    def estim_time(led):

        # === Start
        if led == 0:
            # === Begin time stamp
            led = datetime.datetime.now()
            print('\n\nStart time      : ', led, '\n')

            return led

        # === End of simulation
        else:
            print('\n\nRunning time     : ', datetime.datetime.now() - led)
            print('Main done')


    # === time
    led = 0
    led = estim_time(led)

    # === main
    do_all().run()

    # === time
    led = estim_time(led)
