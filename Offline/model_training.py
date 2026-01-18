import torch
import torch.nn as nn
import torch.nn.functional as F
import random
import numpy as np
import argparse
from torch_geometric.nn import MessagePassing

seed=2022
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

parser = argparse.ArgumentParser(description='Offline Parameters')
parser.add_argument('--dataset','-n',type=str,default='../Dataset',required=False)
parser.add_argument('--dimension','-d',type=int,default=2,required=False)
parser.add_argument('--epoch','-e',type=int,default=1000,required=False)
parser.add_argument('--batch','-b',type=int,default=4096,required=False)
parser.add_argument('--rate','-l',type=float,default=0.01,required=False)

args=parser.parse_args()

path_dataset=str(args.dataset)+'/'
spur_dim=args.dimension
epochs=args.epoch
batch_size=args.batch
learning_rate=args.rate

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class SpanConv(MessagePassing):
    def __init__(self,num_labels,spur_dim=2):
        super().__init__(aggr='add') 
        self.spur_emb = nn.Embedding(num_labels, spur_dim)
    
    def forward(self, x, edge_index):
        spur = F.softplus(self.spur_emb(x))
        span = self.propagate(edge_index, x=spur)
        return spur, span
    
    def message(self, x_j):
        return x_j
    
    def update(self, aggr_out):
        return aggr_out
    
    def forward_spur_only(self,x):
        spur = F.softplus(self.spur_emb(x))
        return spur

def anti_dominance_loss_sampled(z, num_pairs=4096, tau=0.1):
    N, D = z.shape
    idx1 = torch.randint(0, N, (num_pairs,), device=z.device)
    idx2 = torch.randint(0, N, (num_pairs,), device=z.device)

    diff = z[idx1] - z[idx2]
    min_diff = diff.min(dim=-1).values
    soft_dom = torch.sigmoid(min_diff / tau)
    return soft_dom.mean()

data_graph_pyg=torch.load(path_dataset+'data_graph.pt')
num_labels = data_graph_pyg.x.max().item() + 1
print(f"Graph: {data_graph_pyg.num_nodes} nodes, {data_graph_pyg.num_edges} edges")
print(f"Number of classes: {num_labels}")
data_graph_pyg = data_graph_pyg.to(device)

model = SpanConv(num_labels, spur_dim).to(device)
opt = torch.optim.Adam(model.parameters(), lr=learning_rate)

print("Training started...")
for ep in range(1, epochs+1):
    model.train()
    opt.zero_grad()

    spur, span = model(data_graph_pyg.x.squeeze().long(), data_graph_pyg.edge_index)

    loss=anti_dominance_loss_sampled(span,tau=0.5)
    loss.backward()
    opt.step()
    if ep % 100 == 0 or ep == 1:
        print(f"Epoch {ep}/{epochs}, Loss: {loss.item():.6f}")
print("Training finished.")

model.eval()
x=torch.arange(num_labels).long().to(device)
with torch.no_grad():
    spur = model.forward_spur_only(x)

spur_emb = spur.detach().cpu().numpy()
num_vertex, spur_dim = spur_emb.shape
filename = path_dataset + 'spur_emb.txt'
with open(filename, 'w') as f:    
    for i in range(num_vertex):
        line = " ".join(f"{x}" for x in spur_emb[i])
        f.write(line + "\n")

def l1_normalize_rows(z: torch.Tensor, eps: float = 1e-12) -> torch.Tensor:
    row_sum = z.sum(dim=1, keepdim=True).clamp_min(eps)
    return z / row_sum

spur_norm = l1_normalize_rows(spur)
spur_emb = spur_norm.detach().cpu().numpy()
num_vertex, spur_dim = spur_emb.shape
filename = path_dataset + 'spur_emb_norm.txt'
with open(filename, 'w') as f:    
    for i in range(num_vertex):
        line = " ".join(f"{x}" for x in spur_emb[i])
        f.write(line + "\n")
