import torch.nn as nn
H = 64
n_input = 8
n_output = 8
model = nn.Sequential(nn.Linear(n_input, H), nn.ReLU(), nn.BatchNorm1d(H), nn.Dropout(0.2),
                      nn.Linear(H, 2*H), nn.ReLU(), nn.BatchNorm1d(2*H), nn.Dropout(0.3),
                      nn.Linear(2*H, H), nn.ReLU(), nn.BatchNorm1d(H),
                      nn.Linear(H, int(H/2)), nn.ReLU(),
                      nn.Linear(int(H/2), n_output))
