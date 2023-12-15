%% Make Lambda Matrix Function
% Reference : Kyu Ho Kang(2012), Diebold and Li(2006) 
% [Input]

% Yt(tau) = Beta1(t) + Beta2(t)*((1-exp(-lambda(t)*tau))/(lambda(t)*tau)) +
% Beta3(t)*((((1-exp(-lambda(t)*tau))/(lambda(t)*tau)))-exp(-lambda(t)*tau))

% tau : Period Vector
% lambda : Given constant (In K.H Kang(2012) => 0.15)

% [Output]
% Big Lambda Matrix (Lambda)

% Example Useage
%lambda = 0.15;
%tau = [3 6 9 12 15 18 21 24 30 36 48 60 72 84 96 108 120];
% Lambda_M = makeLambda(lambda, tau);

function [Lambda_M] = makeLambda(lambda,tau)

T = length(tau);

Lm = ones(T,3);

for i = 1:T
    for t_para = tau(i)
        Lm(i,2) = (1-exp(-lambda*t_para))/(t_para*lambda);
    end
end

for i = 1:T
    for t_para = tau(i)
        Lm(i,3) = (((1-exp(-lambda*t_para))/(t_para*lambda))-exp(-lambda*t_para));
    end
end

Lambda_M = Lm;

end

% Check the Beta(1,2,3) Loadings
%y1 = Lm(:,1);
%plot(y1)
%hold on
%y2 = Lm(:,2);
%plot(y2);
%y3 = Lm(:,3);
%plot(y3);
%hold off