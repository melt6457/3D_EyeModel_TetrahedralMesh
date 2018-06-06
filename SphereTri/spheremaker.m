function Vert = spheremaker
N = ceil(0.08*[0,14,44,92,184,308,464,704,994,1334,1808,2354]);
Vert = [];
Face = [];

for i = 2:12
    [vMat, fMat] = spheretri(N(i));
    Vert = [Vert;i*vMat];
    Face = [Face;fMat];
end
