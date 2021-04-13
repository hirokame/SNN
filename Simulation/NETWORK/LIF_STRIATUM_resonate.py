#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''Leaky Integrate and Fire neuron Striatum

クラス内変数
R: 膜抵抗
C: 膜容量
tau: 膜電位変化の時定数
uR: 静止電位
thrs: 発火閾値
maxV: 活動電位
v: 膜電位
resv:振動させたあとの膜電位
phase_ipsi: 同側のCortexに対して反応する位相(0~0.9の10段階)
phase_contra: 対側のCortexに対して反応する位相(0~0.9の10段階)
res_int: 閾値を振動させる周期(resonate interval)
res_amp: 閾値を振動させる幅
theta: 振動状態p=[0.3,0.7]を表す角度
res: 振動の大きさ
dop: D1 or D2受容体発現細胞

クラス内メソッド
step: シミュレーションの1ステップ分の操作。method0と1でオイラー法とルンゲクッタの2種類を実装
solve_eular: オイラー法を用いた微分方程式の解法
solve_rk4: ルンゲクッタ法を用いた微分方程式の解法
fv: tau*dV/dt = -V+RIをdVについて解いた式  dV = (-V+RI)/tau*dt

'''
import random
import numpy as np
pi = np.pi

class LIF_striatum:
    def __init__(self):
        self.R = 5  
        self.C = 3
        self.tau = self.R*self.C
        self.uR = -40
        self.rest = -65
        self.thrs = 10
        self.maxV = 90
        self.v = -65
        self.resv = -65
        self.phase_ipsi = random.randint(0, 9)*0.1
        self.phase_contra = random.randint(0, 9)*0.1
        self.res_frq = random.choice([2,3,4])
        self.interval = 1000/self.res_frq
        self.theta = -pi
        self.ref = 10
        self.res = np.sin(self.theta)
        self.dop = random.randint(1,2)
        self.t = 0

    def step(self, dt, I, method=0):
     
        self.t += dt
        dw = (dt/self.interval)*2*pi
        # dw = (dt/self.res_int)*2*pi
        self.theta += dw
        #amp = np.exp(-self.t/1500) # resonantの振幅は1000msで(1/e)倍になる設定
        amp = 1                     # 振幅は減らずに位相だけリセットされるモデル
        self.res = np.cos(self.theta)*amp # １乗
        if self.resv >= self.thrs:
            self.v = self.uR
            self.resv = self.uR
            if self.t > self.interval/2:
                self.theta = 0 ##一定の期間(resonant周期の半分)が経ったあと発火したらresonantをリセット
                self.t = 0
                
        elif self.t < self.ref:
            pass
       
        else: 
            self.resv = self.rest + abs(self.v-self.rest)*self.res # 静止膜電位からずれてる分をresonateした値を静止膜電位に足している
            if method == 1: 
                self.solve_rk4(dt, I)
            elif method == 0: 
                self.solve_euler(dt, I)
            if self.resv >= self.thrs: 
                self.v = self.maxV
                self.resv = self.maxV

    def solve_euler(self, dt, I): #Eular Method
        dv = self.fv(self.v, I) * dt 
        self.dv = dv
        self.v += dv 
        #self.resv = self.rest + abs(self.v-self.rest)*self.res
        self.resv += dv

    def solve_rk4(self, dt, I): #Runge-Kutta Method
        dv1 = self.fv(self.v, I) * dt
        dv2 = self.fv(self.v + dv1 * 0.5, I) * dt
        dv3 = self.fv(self.v + dv2 * 0.5, I) * dt
        dv4 = self.fv(self.v + dv3, I) * dt
        dv = 1/6 * (dv1 + dv2 * 2 + dv3 * 2 + dv4)
        self.v += dv
        #self.resv = self.rest + abs(self.v-self.rest)*self.res
        self.resv += dv

    def fv(self, v, I):
        return (-v+self.R*I)/(self.tau)  ##self.R*self.C = tau