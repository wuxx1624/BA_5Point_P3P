function [im, data3D] = get_patch_ohio(px)

load './data/bigim';

im = I(px.v(1):px.v(2), px.u(1):px.u(2));