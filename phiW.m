% Returns the phase of the DTFT of the window function at angular frequency w

function phase = phiW(w, a, M)

W = sin(M * w / 2) * exp(-1i * w * (M - 1) / 2) * ((a / sin(w / 2)) ...
    - ((1 - a) * ((exp(-1i * pi / M) / sin((w - (2 * pi / M)) / 2)) ...
    + (exp(1i * pi / M) / sin((w + (2 * pi / M)) / 2))) / 2));

phase = angle(W);