#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''Leaky Integrate and Fire neuron GPe

クラス内変数
R: 膜抵抗
C: 膜容量
tau: 膜電位変化の時定数
uR: 静止電位
thrs: 発火閾値
maxV: 活動電位
v: 膜電位
phase: 反応する位相(0~0.9の10段階)

クラス内メソッド
step: シミュレーションの1ステップ分の操作。
solve_eular: オイラー法を用いた微分方程式の解法
solve_rk4: ルンゲクッタ法を用いた微分方程式の解法
fv: tau*dV/dt = -V+RIをdVについて解いた式  dV = (-V+RI)/tau*dt

'''
import random

class LIF:
    def __init__(self):
        self.R = 5  
        self.C = 3
        self.tau = self.R*self.C
        self.uR = -40
        self.thrs = 30
        self.maxV = 90
        self.v = -65
        self.phase = random.randint(0,9)*0.1

    def step(self, I, dt, method=0):  # method==0:Eular Method / method==1:Runge-Kutta Method
        
        if self.v >= self.thrs:
            self.v = self.uR
        else: 
            if method == 1: 
                self.solve_rk4(dt, I)
            elif method == 0: 
                self.solve_euler(dt, I)

            if self.v >= self.thrs: 
                self.v = self.maxV 

    def solve_euler(self, I, dt): #Eular Method
        dv = self.fv(self.v, I) * dt 
        self.v += dv 

    def solve_rk4(self, I, dt): #Runge-Kutta Method
        dv1 = self.fv(self.v, I) * dt
        dv2 = self.fv(self.v + dv1 * 0.5, I) * dt
        dv3 = self.fv(self.v + dv2 * 0.5, I) * dt
        dv4 = self.fv(self.v + dv3, I) * dt
        dv = 1 / 6 * (dv1 + dv2 * 2 + dv3 * 2 + dv4)
        self.v += dv


    def fv(self, v, I):
        return (-v + self.R * I) / (self.tau) ##self.R*self.C = tau

