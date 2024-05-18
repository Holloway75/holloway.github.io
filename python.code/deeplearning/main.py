import os
import logging
import pandas as pd
import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from datetime import datetime
import matplotlib.pyplot as  plt
import seaborn as sns


def masked_mse_loss(y_pred, y_true, mask):
    diff_squared = (y_pred - y_true) ** 2
    masked_diff_squared = diff_squared * mask.float()  # 将布尔掩码转换为浮点数（True -> 1.0, False -> 0.0）

    # 计算有效的元素数量（即掩码中为True的数量）
    num_elements = mask.sum().item()

    # 计算平均误差（注意要避免除以0）
    if num_elements == 0:
        return torch.tensor(0.0, requires_grad=True)  # 如果没有有效的元素，返回0（或任何其他合适的值）
    else:
        mse = masked_diff_squared.sum() / num_elements
        return mse


class NpDataset(Dataset):
    def __init__(self, data):
        """
        Args:
            data (np.ndarray): NumPy数组，包含数据
          """
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        sample = self.data[idx]
        sample = torch.tensor(sample, requires_grad=True, dtype=torch.float32)
        return sample


class Autoencoder(nn.Module):
    def __init__(self, n_traits):
        super().__init__()
        self.n_traits = n_traits
        self.encoder = nn.Sequential(
            nn.Linear(self.n_traits, 200, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(200, 120, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(120,60, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(60, 32, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(32, 4, dtype=torch.float32),
        )
        self.decoder = nn.Sequential(
            nn.Linear(4, 32, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(32, 60, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(60, 120, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(120, 200, dtype=torch.float32),
            nn.ReLU(),
            nn.Linear(200, self.n_traits, dtype=torch.float32)
        )

    def forward(self, x):
        _encoded = self.encoder(x)
        _decoded = self.decoder(_encoded)
        return _encoded, _decoded


def checkpoint_save(filename, _model, _optimizer):
    # 保存训练参数及状态
    check_point = {
        'model_state_dict': _model.state_dict(),
        'optimizer_state_dict': _optimizer.state_dict(),
    }
    torch.save(check_point, filename)


def checkpoint_load(filename, _model, _optimizer):
    _checkpoint = torch.load(filename)

    # 恢复模型的权重和优化器的状态
    _model.load_state_dict(_checkpoint['model_state_dict'])
    _optimizer.load_state_dict(_checkpoint['optimizer_state_dict'])


if __name__ == '__main__':
    os.chdir(r"D:/")
    # 初始训练设置
    epoch_times = 50
    batch_size = 2000
    lr = 0.0001
    torch.set_num_threads(8)
    data = np.load(r"D:\ac10.special.npy", mmap_mode='r')
    dataset = NpDataset(data)

    # 建立模型、优化器、加载器
    model = Autoencoder(n_traits=data.shape[1]).to('cpu')
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    # 加载检查点
    checkpoint_load('checkpoint.50.pth', model, optimizer)

    with open('ac500.samples') as f:
        samples = [i.split("-")[0] for i in f.read().split()]
    data = torch.tensor(data, dtype=torch.float32, requires_grad=True)
    data = model.encoder(data)
    data = data.detach().numpy()
    df = pd.DataFrame(index=samples, columns=[f'v{i}' for i in np.arange(1, 5)], data=data)
    df.to_excel('encoded.xlsx', index=True)
    exit()

    # 配置日志
    logging.basicConfig(filename='training.log', level=logging.INFO)

    loss_fn = nn.MSELoss(reduction='mean')
    logging.info(f'训练开始时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    epoch_losses = []
    for epoch in np.arange(0, epoch_times):
        # 训练
        batch_losses = []
        for batch_samples in dataloader:
            encoded, decoded = model(batch_samples)
            loss = loss_fn(batch_samples, decoded)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            batch_losses.append(loss.detach().numpy())
        epoch_losses.append(np.mean(batch_losses))
        print(f'epoch {epoch}, loss {epoch_losses[-1]}')
        if epoch % 50 == 0 and epoch != 0 :
            checkpoint_save(f'checkpoint.{epoch}.pth', model, optimizer)
    logging.info(f'训练结束时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    for i, loss in enumerate(epoch_losses):
        logging.info(f'epoch {i} 平均损失 {loss}')

    # 保存检查点
    checkpoint_save('check_point.final.pth', model, optimizer)

    # 绘制epoch_loss变动
    fig, ax = plt.subplots()
    sns.lineplot(x=np.linspace(0, 10, epoch_times), y=epoch_losses)
    plt.show()

