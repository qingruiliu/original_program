import tifffile
import numpy as np

def process_tiff(input_file, output_file, start_frame, downsample_rate=3):
    """
    处理 TIFF 序列：
    1. 读取并转化为 grayscale
    2. 从指定帧开始截取
    3. 按指定速率下采样
    4. 保存结果为新的 tiff 序列
    """

    # 1. 读取 tiff 序列
    with tifffile.TiffFile(input_file) as tif:
        images = tif.asarray()

    print(f"原始序列形状: {images.shape}")  # (frames, height, width)

    # 如果是彩色图像 (frames, height, width, channels)，转为灰度
    if images.ndim == 4 and images.shape[-1] in [3, 4]:  
        # 转 grayscale (简单平均法，也可用更复杂的加权)
        images = images.mean(axis=-1).astype(images.dtype)

    # 2. 手动输入的起始帧 (Python 下标从 0 开始)
    if start_frame < 0 or start_frame >= images.shape[0]:
        raise ValueError("起始帧超出范围！")
    images = images[start_frame:]

    print(f"截取后的序列形状: {images.shape}")

    # 3. 下采样 (保留第 0,3,6... 帧)
    images = images[::downsample_rate]

    print(f"下采样后的序列形状: {images.shape}")

    # 4. 保存为新的 tiff 序列
    tifffile.imwrite(output_file, images, dtype=images.dtype)

    print(f"处理完成，保存到: {output_file}")


if __name__ == "__main__":
    # ===== 手动设置 =====
    input_file = "input_stack.tif"       # 输入文件路径
    output_file = "output_downsampled.tif"  # 输出文件路径
    start_frame = int(input("请输入起始帧编号 (从 0 开始): "))  

    process_tiff(input_file, output_file, start_frame, downsample_rate=3)
