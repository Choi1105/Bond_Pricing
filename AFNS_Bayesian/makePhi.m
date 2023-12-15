%  Making Phi matrix
function [Phi, Phi_f, Phi_m] = makePhi(G, GQ)

if rows(G) > 3

    G_ff = G(1:3, 1:3);
    G_fm = G(1:3, 4:end);
    
    Phi_m = G_fm;
    Phi_f = G_ff - GQ;
    Phi = [Phi_f, Phi_m];

else
    
    Phi = G - GQ;

end

end