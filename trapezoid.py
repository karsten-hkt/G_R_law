import numpy as np
import matplotlib.pyplot as plt

# 定义梯形源时函数的参数
start_time = 0.0      # 起始时间
rise_time = 2.0       # 上升时间
top_time = 4.0        # 平顶时间
fall_time = 3.0       # 下降时间
amplitude = 1.0       # 峰值幅度

# 生成时间序列
total_time = start_time + rise_time + top_time + fall_time
time = np.linspace(0, total_time, 1000)

# 计算梯形源时函数
source_time_function = np.zeros_like(time)
source_time_function[(time >= start_time) & (time < start_time + rise_time)] = (time[(time >= start_time) & (time < start_time + rise_time)] - start_time) / rise_time * amplitude
source_time_function[(time >= start_time + rise_time) & (time < start_time + rise_time + top_time)] = amplitude
source_time_function[(time >= start_time + rise_time + top_time) & (time <= total_time)] = amplitude - (time[(time >= start_time + rise_time + top_time) & (time <= total_time)] - (start_time + rise_time + top_time)) / fall_time * amplitude

# 绘制梯形源时函数
plt.figure(figsize=(8, 4))
plt.plot(time, source_time_function, label='Trapezoid Source Time Function')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Trapezoid Source Time Function')
plt.grid(True)
plt.legend()
plt.show()
