import os
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from rastermap import utils, Rastermap

# 生成一些示例数据
np.random.seed(42)
num_neurons = 100
num_timepoints = 1000
data = np.random.rand(num_neurons, num_timepoints) # 100个神经元，每个神经元有1000个时间点
#重新随机排列神经元和时间点的顺序
data = data[:, np.random.permutation(num_timepoints)]

# 初始化Rastermap对象
rastermap = Rastermap(n_components=2)  # 使用正确的初始化参数

# 拟合数据并进行降维
embedding = rastermap.fit_transform(data)

# 可视化结果
plt.figure(figsize=(10, 6))
plt.scatter(embedding[:, 0], embedding[:, 1], c=np.arange(num_neurons), cmap='viridis')
plt.colorbar(label='Neuron Index')
plt.xlabel('Component 1')
plt.ylabel('Component 2')
plt.title('Rastermap Embedding')
plt.show()

#https://github.com/MouseLand/course-materials/blob/main/behavior_encoding/tutorial_solutions.ipynb
#https://github.com/MouseLand/course-materials/blob/main/behavior_encoding/tutorial.ipynb