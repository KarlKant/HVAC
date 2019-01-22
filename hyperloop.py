#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 10:02:43 2018

@author: Xiebingchen
"""
import math,sys,random
import matplotlib.pyplot as plt
def progress_bar(num,tot):
    rate = float(num)/tot
    ratenum = int(100*rate)
    r = '\r[{}{}]{}%'.format('*'*ratenum,' '*(100-ratenum),ratenum)
    sys.stdout.write(r)
    sys.stdout.flush()
    
#基本常量
Pai=3.1415926#圆周率
Rou=1.293#常温常压空气密度
Sita=5.67e-8#玻尔兹曼常量
R=287#R常量
C_water=4200
Cv=2.5*R#定容比热容
Cp=3.5*R#定压比热容
kao=Cp/Cv

#几何对象
class cylinder:
    '套筒圆筒基类'
    def __init__(self,r,h,T,Yipu):
        self.r=r
        self.h=h
        self.T=T+273.15#表面温度
        self.Yipu=Yipu
        self.V=Pai*r*r*h
        self.S=Pai*r*2*h
        self.names={'r':'半径',
                    'h':'高度',
                    'T':'温度',
                    'Yipu':'发射率',
                    'V':'体积',
                    'S':'侧面积'}
        self.dimension={'r':'m',
                        'h':'m',
                        'T':'K',
                        'Yipu':'',
                        'V':'m^3',
                        'S':'m^3'}
#换热对象
class Heatchange:
    '换热基类'
    def __init__(self,d,T1,T2,lanb,alpha):
        #d平板厚度，T1对流侧温度，T2辐射侧温度，lanb热导率，alpha热扩散系数
        self.d=d#平板厚度
        self.T1=T1#一侧温度
        self.T2=T2#另一侧温度
        self.lanb=lanb#热导率
        self.alpha=alpha#热扩散系数
        self.names={}#
        self.dimension={}#
        
    '''边界一侧是对流，另一侧是辐射'''
    def unstable(self,d_step,time_step,time,T0,h,T_,S1,S2,Yipu1,Yipu2):
        #d_step空间步长，time_step时间步长，time总时间，T0对流温度，h对流系数
        #T_辐射温度，S1外界面积，S2平板面积，Yipu1外界辐射度，Yipu2平板辐射度
        temps=[]
        temp=[]
        temp0=[]
        Fo=self.alpha*time_step/(d_step*d_step)
        if Fo<0.5:
            Bi=h*d_step/self.lanb
            for x in range(0,int(self.d/d_step)+1):
                temp.append(self.T1)
            temps.append(temp.copy())
            for t in range(1,int(time/time_step)+1):
                temp_last=temp.copy()
                for y in range(0,int(self.d/d_step)+1):
                    if y==0:
                        temp[y]=temp_last[y]*(1-2*Fo*Bi-2*Fo)+2*Fo*temp_last[y+1]+2*Fo*Bi*T0
                    elif y==int(self.d/d_step):
                        radiation=sleeve_radiation(T_,temp_last[y],S1,S2,Yipu1,Yipu2)#平板为S2,外部为S1
                        temp[y]=(temp_last[y-1]-temp_last[y])*Fo+temp_last[y]+2*time_step*self.alpha*radiation/(self.lanb*d_step*S1)
                    else:
                        temp1=temp_last[y+1]
                        temp2=temp_last[y-1]
                        temp3=temp_last[y]
                        temp[y]=self.alpha * time_step * (temp1 + temp2) / (d_step * d_step) + (1 - 2 * self.alpha * time_step / (d_step * d_step)) * temp3
                if t*time_step==int(t*time_step):
                    temps.append(temp.copy())
                if t*time_step/100==int(t*time_step/100):
                    progress_bar(t*(self.d/d_step+1)+y+1,time/time_step*(self.d/d_step+1))
        else:
            return('Fol_error')
        for i in range(0,int(time)+1):
            for j in range(0,int(self.d/d_step)+1):
                temps[i][j]=round(temps[i][j]-273.15,2)
                if j==0:
                    temp0.append(temps[i][j])
        return {'temps':temps,
                'temp0':temp0}
        
#辐射换热函数
def sleeve_radiation(T1,T2,S1,S2,Yipu1,Yipu2):
        E1=Sita*math.pow(T1,4)
        E2=Sita*math.pow(T2,4)
        return (E1-E2)/(1/(Yipu1*S1)+1/(Yipu2*S2)-1/S2)

        
#空调换热函数
'''基本思路：采用变风量送风系统
将冷负荷分成预估部分和偏差部分两个部分
预估冷负荷由计算获取，负责决定不同时间的基础送风量
分为辐射换热和内部人员、灯管三种负荷
其中车体的辐射换热还要考虑车体的导热和内部的对流，当内部温度控制在26+-0.5摄氏度时，理论与实际偏差不大
偏差部分根据测得实际的室内温度获得，对送风量进行修正
'''
#负荷送风量：送风温度,实际温度,负荷-->输出送风量(单位时间)
def loadair(T_inlet,T_want,load):
    m_air=load/((T_want-T_inlet)*Cp)#一秒送气质量
    return m_air/Rou
    

#修正送风量：送风温度，室内温度，设计温度，空间总体积-->输出送风量（总 修正偏差所需）
def coolair(T_inlet,T_real,T_want,V_total):
    return (T_real-T_want)*V_total/(T_want-T_inlet)
    
#实际送风量
class fengji:
    '风机基类'
    def __init__(self,Q,T):#送风半径，最大送风速度，送风温度//最大送风量，送风温度
        self.Q=Q
        self.T=T+273.15
    
#实际温度计算
class realT:
    '室内温度基类'
    def __init__(self,T0,V_total):
        self.T_real=T0
        self.V_total=V_total
        self.T=[T0,]
    def realTemp(self,load,Q_air,T_air):
        P_heat=load+(T_air-self.T_real)*Cp*Q_air*Rou
        dT=P_heat/(Cp*self.V_total*Rou)
        self.T_real+=dT
        self.T.append(self.T_real)
        return 0
    def tempList(self):
        T_=[]
        for i in self.T:
            i-=273.15
            T_.append(i)
        return T_

#换气对象
class airChange:
    '换气基类'
    def __init__(self,changeRate,totalV):
        self.changeRate=changeRate
        self.totalV=totalV
        self.V=changeRate*totalV/3600
        self.m=changeRate*totalV*Rou/3600
        self.names={'Rou':'密度','changeRate':'换气次数','totalV':'空间总体积','V':'换气体积','m':'换气质量'}
        self.dimension={'Rou':'千克每立方米','changeRate':'次每小时','totalV':'立方米','V':'立方米每秒','m':'千克每秒'}

#送水回水函数
def water(bool_hc,Q_air,T_air):
    if bool_hc=='cool':
        waterT1=5+273.15
        waterT2=12+273.15
        airT1=T_air
        airT2=cold_air+273.15
        veryT1=125
        veryT2=300
    if bool_hc=='heat':
        waterT1=45+273.15
        waterT2=40+273.15
        airT1=T_air
        airT2=heat_air+273.15
        veryT1=857
        veryT2=300
    Q_water=-Cp*Rou*Q_air*(airT2-airT1)/(C_water*(waterT2-waterT1)*YiTa_water)
    Q_very=-Q_air*(airT2-airT1)/(veryT2-veryT1)/Rou/YiTa_water/YiTa_water
    return [Q_water,Q_very]
        





#输入变量
'''设计输入变量'''
h_air=5
'''隧道输入参数'''
r_tube=1.83#管道半径
YiPu_tube=0.85#管道辐射
'''车厢输入参数'''
d_wall=0.2#车壁厚度
r_pod=1.25#车厢半径
l_pod=25#车厢长度
YiPu_pod=0.55#车厢辐射度
alpha_pod=237/(2700*0.88*103)#热扩散系数
lanbu_pod=237#热导率
'''非稳态计算输入'''
time_step=0.01
long_step=0.01
time_long=3600
'''其它参数'''
YiTa_water=0.7

#可修改输入变量
'''设计输入变量'''
changeRate=2#换气次数
temp_atm=26#大气气温
'''隧道输入参数'''
T_tube=50#管道内壁温度
'''车厢输入参数'''
temp_wall=temp_atm
temp_in=26#车厢内部温度设定
'''送风输入'''
cold_air=20
heat_air=28
Q_max=10
'''随机波动'''
randoms=2000


print('===预计算===')
#定义对象
c_tube=cylinder(r_tube,l_pod,T_tube,YiPu_tube)#半径1.83，长度25，内壁温度T_tube，辐射度YiTa_tube
c_pod=cylinder(r_pod,l_pod,temp_wall,YiPu_pod)#半径1.25，长度25，稳定温度26，辐射度YiTa_pod

heatC=Heatchange(d_wall,c_pod.T,c_pod.T,lanbu_pod,alpha_pod)#换热对象
#d平板厚度，T1对流侧温度K，T2辐射侧温度K，lanb热导率，alpha热扩散系数
air=airChange(changeRate,c_pod.V)#换气对象
T_pod=heatC.unstable(long_step,time_step,time_long,temp_in+273.15,h_air,c_tube.T,c_pod.S,c_tube.S,c_pod.Yipu,c_tube.Yipu)#非稳态计算函数
#d_step空间步长，time_step时间步长，time总时间，T0对流温度，h对流系数
#T_辐射温度，S1外界面积，S2平板面积，Yipu1外界辐射度，Yipu2平板辐射度



print(' ')
print('===输入参数===')
print('换气次数',changeRate,'次每小时')
print('车厢内壁对流换热系数',h_air,'')
print('管道内壁温度',T_tube,'摄氏度')
print('管道内壁辐射度',YiPu_tube,'')
print('车体初始温度',temp_wall,'摄氏度')
print('车体表面辐射度',YiPu_pod,'')
print('热扩散系数',alpha_pod,'')
print('热导率',lanbu_pod,'')
print('迭代时间步长',time_step,'秒')
print('迭代空间步长',long_step,'米')
print('总运行时间',time_long,'秒')
print(air.names['totalV'],round(air.totalV,2),air.dimension['totalV'])



'''实际温度计算'''
realtemp=realT(temp_in+273.15,c_pod.V)
flag=0
Q_extra=0
inlet_air_Q=[]
inlet_air_temp=[]
P_convections=[]
P_loads=[]
P_realLoads=[]
fresh_airs=[]
unfresh_airs=[]
waterQc=[]
waterQh=[]
veryQc=[]
veryQh=[]
for t in range(0,time_long+1):
    fresh_airs.append(air.V)
    T_real=realtemp.T_real
    T_want=temp_in+273.15
    #辐射负荷计算，实际上是对流计算，利用前面算出来的列表
    T_convection1=T_pod['temp0'][t]+273.15
    T_convection2=T_real
    P_convection=h_air*(T_convection1-T_convection2)*c_pod.S
    P_convections.append(P_convection)
    #人员负荷计算，一个单变量函数
    P_people=2000*(0.75+0.5*t/(time_long+1))
    #灯光负荷，无变量函数，不算函数，哈哈哈哈哈
    P_light=800#每平米11瓦特
    P_load=P_light+P_convection+P_people
    P_loads.append(P_load)
    #送风量---根据前面两个送风量的计算函数来，送风量需要和风机结合起来
    if flag==0:
        if P_load>0:
            air_inletemp=cold_air+273.15
        else:
            air_inletemp=heat_air+273.15
    if (T_real-T_want)>0.5:
        Q_extra=coolair(air_inletemp,T_real,T_want,c_pod.V)
        flag=1
        Q=Q_max
    elif (T_real-T_want)<(-0.5):   
        Q_extra=coolair(air_inletemp,T_real,T_want,c_pod.V)
        flag=1
        Q=Q_max
    if Q_extra<0:
        flag=2
        if air_inletemp==cold_air+273.15:
            air_inletemp=heat_air+273.15
        else:
            air_inletemp=cold_air+273.15
        if (T_real-temp_in)>0.5:
            Q_extra=coolair(air_inletemp,T_real,T_want,c_pod.V)
        elif (T_real-temp_in)<(-0.5):   
            Q_extra=coolair(air_inletemp,T_real,T_want,c_pod.V)
    Q_main=loadair(air_inletemp,T_want,P_load)
    if flag==1:
        Q_real=Q_max
        Q_extra-=(Q_real-Q_main)
        if Q_extra<0:
            flag=0
            Q_extra=0
    elif flag==2:
        Q_real=Q_max
        Q_extra-=(Q_real+Q_main)
        if Q_extra<0:
            flag=0
            Q_extra=0
    else:
        if Q_main>Q_max:
            Q_real=Q_max
        else:
            Q_real=Q_main
    ###最好可以生成一个List，对风机情况进行记录
    #print('======',len(realtemp.tempList()))
    #print(Q_real,air_inletemp)
    inlet_air_Q.append(Q_real)
    if air_inletemp==cold_air+273.15:
        bool_hc='cool'
        waterQ=water(bool_hc,Q_real,T_real)
        waterQc.append(waterQ[0]*1000)#千克换算立方厘米
        veryQc.append(waterQ[1])
        waterQh.append(0)
        veryQh.append(0)
    else:
        bool_hc='heat'
        waterQ=water(bool_hc,Q_real,T_real)
        waterQh.append(waterQ[0]*1000)#千克换算立方厘米
        veryQh.append(waterQ[1])
        waterQc.append(0)
        veryQc.append(0)
    if (Q_real-air.V)>0:
        unfresh_airs.append(Q_real-air.V)
    else:
        unfresh_airs.append(0)
    inlet_air_temp.append(air_inletemp-273.15)
    P_loadReal=P_load+random.randint(-randoms,randoms)
    P_realLoads.append(P_loadReal)
    realtemp.realTemp(P_loadReal,Q_real,air_inletemp)
tempList=realtemp.tempList()

realoads_av=[]
t_load=0 

for t in range(0,time_long+1):
    averagel=60
    if t/averagel==int(t/averagel):
        if t!=0:
            av_load=t_load/averagel
            for i in range(0,averagel):
                realoads_av.append(av_load)
        t_load=0
    t_load+=P_realLoads[t]
    
    
print('===计算结果===')
plt.figure(1)
plt.plot(T_pod['temp0'])
plt.title('Interior wall temperature')
plt.xlabel('Time/s')
plt.ylabel('Temperature/°C')    
plt.figure(2)
plt.plot(tempList)
plt.title('Inside air temperature')
plt.xlabel('Time/s')
plt.ylabel('Temperature/°C')
plt.figure(3)
plt.plot(inlet_air_Q)
plt.title('Fan flow')
plt.xlabel('Time/s')
plt.ylabel('Q/(m³/s)')
plt.figure(4)
plt.plot(inlet_air_temp)
plt.title('Supply air temperature')
plt.xlabel('Time/s')
plt.ylabel('Temperature/°C')
plt.figure(5)
plt.plot(P_convections)
plt.title('Load of heat transfer')
plt.xlabel('Time/s')
plt.ylabel('Load/W')
plt.figure(6)
plt.plot(realoads_av)
plt.plot(P_loads)
plt.title('Total Load (blue:real;orange:estimated) ')
plt.xlabel('Time/s')
plt.ylabel('Load/W')
plt.figure(7)
plt.plot(fresh_airs)
plt.title('Valve flow of freshair')
plt.xlabel('Time/s')
plt.ylabel('Q/(m³/s)')
plt.figure(8)
plt.plot(fresh_airs)
plt.plot(unfresh_airs)
plt.title('Valve flow (orange:2,blue:1&3)')
plt.xlabel('Time/s')
plt.ylabel('Q/(m³/s)')
plt.figure(9)
plt.plot(waterQc)
plt.plot(waterQh)
plt.title('Water flow')
plt.xlabel('Time/s')
plt.ylabel('Q/(cm³/s)')
plt.figure(10)
plt.plot(veryQc)
plt.plot(veryQh)
plt.title('Cold/Heat source flow')
plt.xlabel('Time/s')
plt.ylabel('Q/(kg/s)')


'''
声速：347.252
运动粘度 161.5
导热系数 0.026
导热系数 0.026
普朗特数 0.719
                                
'''
    
    
  

    
