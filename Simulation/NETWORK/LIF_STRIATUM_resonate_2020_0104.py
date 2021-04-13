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

x,y:Resonateの状態空間上の軸。yがVに対応。
r: 曲座標平面における原点(0,-65)からの距離

'''
import random
import numpy as np
pi = np.pi

class LIF_striatum:
    def __init__(self):
        self.R = 5  
        self.C = 3
        self.uR = 25
        self.rest = 0
        self.thrs = 35
        self.maxV = 155
        self.v = 0
        self.phase_ipsi = random.randint(0, 9)*0.1
        self.phase_contra = random.randint(0, 9)*0.1
        self.res_frq = random.choice([2,3,4])
        self.interval = 1000/self.res_frq
        self.theta = -0.5*pi
        self.ref = 10
        self.dop = random.randint(1,2)
        self.t = 0
        self.tau = 0
        self.x = 0
        self.y = 0
        self.r = 0
        self.count = 0

    def step(self, dt, I, method=0):
     
        self.t += dt
        self.tau += dt
        dw = (dt/self.interval)*2*pi
        # dw = (dt/self.res_int)*2*pi
        self.theta += dw
 
        #amp = np.exp(-self.t/1500) # resonantの振幅は1000msで(1/e)倍になる設定
        amp = 1                    # resonantの減衰はおこらない設定
        self.x = self.r*np.cos(self.theta)*amp
        self.y = self.r*np.sin(self.theta)*amp + self.rest #回転軌道上を1Step動かす
        
        self.v = self.y ## self.vに1Step回転させた後の極座標平面上のyの値を対応させる
        if self.v >= self.thrs: #1Step動かして閾値を超えているかを判定
            self.v = self.uR
            self.tau = 0
            if self.t > self.interval/2: 
                self.theta = 0 ##一定の期間(resonant周期の半分)が経ったあと発火したらresonantをリセット
                self.t = 0
                
        elif self.tau < self.ref:
            pass
            
        else: 
            if method == 1:
                self.solve_rk4(dt, I)
                self.y = self.v
                
                if self.count == 0:
                    self.thrs = 35
                else:
                    self.thrs = 35 - 5*self.count  ##繰り返しに応じて閾値が下がっていく。
                
                self.theta = np.arctan2(self.y, self.x) #self.vが変更された後の曲座標平面でのthetaを計算
                if self.theta < 0:
                    self.theta += 2*pi #-pi~piの表記を0~2*piに変換
                self.r = np.sqrt(self.x**2+self.y**2)
                
            elif method == 0: 
                self.solve_euler(dt, I)
                self.y = self.v
                self.theta = np.arctan2(self.y, self.x) #self.vが変更された後の曲座標平面でのthetaを計算
                if self.theta < 0:
                    self.theta += 2*pi #-pi~piの表記を0~2*piに変換
                self.r = np.sqrt(self.x**2+self.y**2)
                
            if self.v >= self.thrs: 
                self.v = self.maxV
                self.tau = 0
                if 1.8*pi<self.theta<2.2*pi:
                    self.count+=1
                    self.count=min(3,self.count)  ##3回くらいまでは繰り返しで閾値が下がる
                else:
                    self.count=0
                    
                if self.t > self.interval/3:
                    self.theta = 0
                    self.t = 0

    def solve_euler(self, dt, I): #Eular Method
        dv = self.fv(self.v, I) * dt 
        self.dv = dv
        self.v += dv 

    def solve_rk4(self, dt, I): #Runge-Kutta Method
        dv1 = self.fv(self.v, I) * dt
        dv2 = self.fv(self.v + dv1 * 0.5, I) * dt
        dv3 = self.fv(self.v + dv2 * 0.5, I) * dt
        dv4 = self.fv(self.v + dv3, I) * dt
        dv = 1/6 * (dv1 + dv2 * 2 + dv3 * 2 + dv4)
        self.v += dv
        

    def fv(self, v, I):
        return (-v + self.R*I) / (self.R*self.C)  ##self.R*self.C = tau