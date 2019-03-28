import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os.path

# 0 - Формула левых прямоугольников
# 1 - Формула правых прямоугольников
# 2 - Формула средних прямоугольников
# 3 - Формула трапеций
# 4 - Формула Симпсона
# 5 - "Правило трех восьмых"

# Treads,Id,Steps,Time,Integral,AbsErr,RungeErr

if __name__ == "__main__":
    names = [
        "Формула левых прямоугольников",
        "Формула правых прямоугольников",
        "Формула средних прямоугольников",
        "Формула трапеций",
        "Формула Симпсона",
        "Формула трех восьмых",
        "Метод Монте-Карло"]


    if os.path.exists('out1.csv'):
        if not os.path.exists("img"):
            os.mkdir("img")
        raw_data = pd.read_csv('out1.csv')

        data_time = raw_data[['Threads', 'Id', 'Steps', 'Time']].groupby(
            ['Threads', 'Id', 'Steps'], as_index=False).mean()

        max_y_time = max(data_time['Time'])
        max_y_sp = 2
        max_y_ep = 2
        for Id in sorted(set(data_time['Id'])):
            sub_df = data_time.loc[data_time.Id == Id]
            fig = plt.figure(figsize=(8, 6), dpi=100)
            fig.suptitle(names[Id], fontsize=16)
            rt_subplt = fig.add_subplot(211)
            plt.ylim(0, max_y_time)
            plt.grid(ls=':')
            Sp_subplt = fig.add_subplot(223)
            plt.ylim(0, max_y_sp)
            plt.grid(ls=':')
            Ep_subplt = fig.add_subplot(224)
            plt.ylim(0, max_y_ep)
            plt.grid(ls=':')
            rt_subplt.set_title('Время выполнения')
            rt_subplt.set_xlabel('Число потоков')
            rt_subplt.set_ylabel('Время, сек')
            Sp_subplt.set_title('Ускорение')
            Sp_subplt.set_xlabel('Число потоков')
            Sp_subplt.set_ylabel("$S_{p}$")
            Ep_subplt.set_title('Эффективность')
            Ep_subplt.set_xlabel('Число потоков')
            Ep_subplt.set_ylabel("$E_{p}$")

            for Step in sorted(set(sub_df['Steps'])):
                sub_step = sub_df.loc[sub_df.Steps == Step]
                one_thread_t = float(sub_step[sub_step.Threads == 1]['Time'])
                speedup = one_thread_t / np.array(sub_step['Time'])
                efficiency = speedup / np.array(sub_step['Threads'])

                rt_subplt.plot(sub_step['Threads'], sub_step['Time'],
                            marker=".", label=f"{Step}")
                Sp_subplt.plot(sub_step['Threads'], speedup, marker=".",
                            label=f"{Step}")
                Ep_subplt.plot(sub_step['Threads'], efficiency, marker=".",
                            label=f"{Step}")

            rt_subplt.legend()
            fig.subplots_adjust(wspace=0.5, hspace=0.5)
            fig.savefig("img/1d-" + str(Id) + ' - ' + names[Id] + '.png')

    if os.path.exists('out2.csv'):
        if not os.path.exists("img"):
            os.mkdir("img")
        raw_data = pd.read_csv('out2.csv')

        data_time = raw_data[['Threads', 'Id', 'Steps', 'Time']].groupby(
            ['Threads', 'Id', 'Steps'], as_index=False).mean()

        max_y_time = max(data_time['Time'])
        max_y_sp = 2
        max_y_ep = 2
        for Id in sorted(set(data_time['Id'])):
            sub_df = data_time.loc[data_time.Id == Id]
            fig = plt.figure(figsize=(8, 6), dpi=100)
            fig.suptitle(names[Id], fontsize=16)
            rt_subplt = fig.add_subplot(211)
            plt.ylim(0, max_y_time)
            plt.grid(ls=':')
            Sp_subplt = fig.add_subplot(223)
            plt.ylim(0, max_y_sp)
            plt.grid(ls=':')
            Ep_subplt = fig.add_subplot(224)
            plt.ylim(0, max_y_ep)
            plt.grid(ls=':')
            rt_subplt.set_title('Время выполнения')
            rt_subplt.set_xlabel('Число потоков')
            rt_subplt.set_ylabel('Время, сек')
            Sp_subplt.set_title('Ускорение')
            Sp_subplt.set_xlabel('Число потоков')
            Sp_subplt.set_ylabel("$S_{p}$")
            Ep_subplt.set_title('Эффективность')
            Ep_subplt.set_xlabel('Число потоков')
            Ep_subplt.set_ylabel("$E_{p}$")

            for Step in sorted(set(sub_df['Steps'])):
                sub_step = sub_df.loc[sub_df.Steps == Step]
                one_thread_t = float(sub_step[sub_step.Threads == 1]['Time'])
                speedup = one_thread_t / np.array(sub_step['Time'])
                efficiency = speedup / np.array(sub_step['Threads'])

                rt_subplt.plot(sub_step['Threads'], sub_step['Time'],
                            marker=".", label=f"{Step}")
                Sp_subplt.plot(sub_step['Threads'], speedup, marker=".",
                            label=f"{Step}")
                Ep_subplt.plot(sub_step['Threads'], efficiency, marker=".",
                            label=f"{Step}")

            rt_subplt.legend()
            fig.subplots_adjust(wspace=0.5, hspace=0.5)
            fig.savefig("img/2d-" + str(Id) + ' - ' + names[Id] + '.png')
