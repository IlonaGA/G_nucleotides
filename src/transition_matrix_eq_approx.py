#!/usr/bin/env python
# coding: utf-8

# In[9]:


import torch
import numpy as np


# In[16]:


def g(x):
    """
    replace with real differntial version
    """
    mask = torch.abs(x) < 1e-6
    x1 = torch.where(mask, x, 0.0)
    x2 = torch.where(mask, 1.0, x)
    
    return torch.where(mask, 1 / (1 - x1), 2 * x2 / (1 - torch.exp(-2 * x2)))

def transition_matrix(v, w, u, Neff):
    """
    v : [2,4] tensor
    w : [4,4] tensor
    u : [4,4] tensor
    Neff: scalar tensor (effective population size, needs to be fitted like v and w)
    it needs a working g function in scope
    """

    def s(bd, ac):
        """
        Neff * (e(bd) - e(ac)) (both numbers in {0,...,15})
        """
        a, c = ac // 4, ac % 4
        b, d = bd // 4, bd % 4

        s_bd = v[0, b] + v[1, d] + w[b, d]
        s_ac = v[0, a] + v[1, c] + w[a, c]

        return torch.tensor([s_bd - s_ac])

    def q(bc, ac):

        if s(bc, ac) < 0:
            a = ac // 4
            b = bc // 4
            return -u[a, b] / (s(bc, ac) / Neff)
        else:
            return torch.zeros([])

    def rate(ac, bd):
        a, c = ac // 4, ac % 4
        b, d = bd // 4, bd % 4

        # Single mutation
        if a == b and c != d:
            return u[c, d] * g(s(bd, ac))

        # Other single mutation
        if a != b and c == d:
            return u[a, b] * g(s(bd, ac))

        # Double mutation
        return (
            q(b * 4 + c, a * 4 + c) * u[c, d] + q(a * 4 + d, a * 4 + c) * u[a, b]
        ) * g(s(bd, ac))

    return torch.stack([rate(a, b) for a in range(16) for b in range(16)]).reshape(
        [16, 16]
    )


# # Test

# In[17]:


v = torch.rand([2, 4], dtype=torch.float64, requires_grad=True)
w = torch.rand([4, 4], dtype=torch.float64, requires_grad=True)
mutation_rate = 1e-6
u = np.full((4, 4), mutation_rate)
for i in range(4):
    u[i, i] = 1 - 3 * mutation_rate
u = torch.from_numpy(u)
Neff = torch.full([], 1000.0, requires_grad=True, dtype=torch.float64)


# In[18]:


transition_matrix(v, w, u, Neff)
