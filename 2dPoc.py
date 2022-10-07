import torch


def apply_rotation(pos, rot):
    x, y = pos[0], pos[1]
    c, s = rot.cos(), rot.sin()
    # TODO: check if right
    xx = x * c + y * s
    yy = x * s - y * c
    ret = torch.stack([xx, yy])
    assert ret.shape == pos.shape
    return ret


class Bead:
    def __init__(self, pos):
        self.pos = pos


class Chain:
    def __init__(self, beads=None):
        self.beads = beads or []


class Block:
    def __init__(self, pos, rot, glued_beads=None):
        self.pos = pos
        self.rot = rot
        self.glued_beads = glued_beads or []  # list of (Bead, relative_position)

    def update_glued_beads(self):
        for b, rel_pos in self.glued_beads:
            b.pos = apply_rotation(rel_pos, self.rot) + self.pos


variables = []
blocks = []
beads = []
chains = []

link_length = 20
dist_k = 1.0 / (link_length**2)
angle_k = 2.0

import math


def calculate_loss():
    for block in blocks:
        block.update_glued_beads()

    losses = []
    for chain in chains:
        for i in range(len(chain.beads) - 1):
            a = chain.beads[i].pos
            b = chain.beads[i + 1].pos
            diff = a - b
            # TODO: proper way to get squared norm
            dist = diff.norm()
            prob = dist - link_length
            losses.append(prob * prob * dist_k)

        for i in range(len(chain.beads) - 2):
            a = chain.beads[i].pos
            b = chain.beads[i + 1].pos
            c = chain.beads[i + 2].pos
            ab = a - b
            cb = c - b
            ab_angle = torch.atan2(ab[1], ab[0])
            cb_angle = torch.atan2(cb[1], cb[0])
            diff = ab_angle - cb_angle
            abs_diff = diff.abs()
            from_straight = math.pi - abs_diff
            losses.append(from_straight * from_straight * angle_k)
    return torch.stack(losses).sum()


center = torch.tensor([250.0, 250.0])

import random


def make_block():
    pos = torch.autograd.Variable(torch.normal(center, 20.0), requires_grad=True)
    rot = torch.autograd.Variable(torch.tensor(3.0), requires_grad=True)
    block = Block(pos, rot)
    variables.append(pos)
    variables.append(rot)
    blocks.append(block)
    return block


def make_glued_bead():
    pos = torch.autograd.Variable(torch.normal(center, 20.0))
    bead = Bead(pos)
    beads.append(bead)
    return bead


def make_floating_bead():
    pos = torch.autograd.Variable(torch.normal(center, 50.0), requires_grad=True)
    bead = Bead(pos)
    beads.append(bead)
    variables.append(pos)
    return bead


def make_chain(l=10):
    chain = Chain()
    chains.append(chain)
    for _ in range(l):
        chain.beads.append(make_floating_bead())
    return chain


make_chain(20)
# make_chain(15)
make_chain(5)

block1 = make_block()
block2 = make_block()

bead1 = make_glued_bead()
block1.glued_beads.append((bead1, torch.tensor([30, 12])))
chains[0].beads.append(bead1)

bead1 = make_glued_bead()
block1.glued_beads.append((bead1, torch.tensor([20, 12])))
chains[0].beads.append(bead1)

bead1 = make_glued_bead()
block2.glued_beads.append((bead1, torch.tensor([-20, 12])))
chains[0].beads.insert(0, bead1)

bead1 = make_glued_bead()
block2.glued_beads.append((bead1, torch.tensor([-20, 22])))
chains[0].beads.insert(1, bead1)


bead1 = make_glued_bead()
block2.glued_beads.append((bead1, torch.tensor([30, 12])))
chains[1].beads.append(bead1)

bead1 = make_glued_bead()
block2.glued_beads.append((bead1, torch.tensor([20, 12])))
chains[1].beads.append(bead1)

bead1 = make_glued_bead()
block1.glued_beads.append((bead1, torch.tensor([-20, 12])))
chains[1].beads.insert(0, bead1)

bead1 = make_glued_bead()
block1.glued_beads.append((bead1, torch.tensor([-20, 22])))
chains[1].beads.insert(1, bead1)


import pygame

pygame.init()

# create

clock = pygame.time.Clock()
screen = pygame.display.set_mode([500, 500])


optimizer = torch.optim.Adam(variables, lr=1)
# optimizer = torch.optim.SGD(variables, lr=0.1, momentum=0.9)


def tensor_to_tuple(a):
    return tuple(map(int, a.detach().numpy()))


t = 0
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill((255, 255, 255))

    for block in blocks:
        w = 60
        h = 40
        points = [(w / 2, h / 2), (-w / 2, h / 2), (-w / 2, -h / 2), (w / 2, -h / 2)]
        points = [
            apply_rotation(torch.tensor(p), block.rot) + block.pos for p in points
        ]
        print(points[0])
        points = [tuple(map(int, p.detach().numpy())) for p in points]
        print(points)
        pygame.draw.polygon(screen, (0, 128, 0), points)
        pygame.draw.circle(screen, (255, 0, 255), tensor_to_tuple(block.pos), 3)

    for chain in chains:

        for i in range(len(chain.beads) - 1):
            a = chain.beads[i].pos
            b = chain.beads[i + 1].pos
            pygame.draw.line(
                screen, (0, 0, 255), tensor_to_tuple(a), tensor_to_tuple(b), 2
            )

        for bead in chain.beads:
            pygame.draw.circle(screen, (0, 0, 255), tensor_to_tuple(bead.pos), 3)

    pygame.display.flip()

    t += 1
    for _ in range(10):
        optimizer.zero_grad()
        loss = calculate_loss()
        loss.backward()
        optimizer.step()
        # print(block1.pos.grad)

    if t % 30 == 0:
        for v in beads:
            v.pos.data = torch.normal(v.pos.data, 10)
        print("asdf")

    clock.tick(30)

pygame.quit()
