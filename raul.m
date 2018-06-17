

g = @(qz, omm, om) (1./sqrt(omm.*(1 + qz).^3 + (1 - omm).*(1 + qz).^(3 + 3.*om)));
distluz = @(omm, om, z) (1 + z).*(integral(@(qz) g(qz, omm, om), 0, z));

z = [0 1 5];
dist = [0 0 0];
for i = 1:length(z)
        dist(i) = distluz(0.3, -1, z(i));
end
%dist = distluz(0.3, -1, z);

plot(z, dist);

pedido = distluz(0.3, -1, 0.5);
disp(fprintf('El valor pedido es: %f', pedido)) ;

load('Union2.mat');

om = -1;
omm = 0.3;
numM = zeros(1, 557);
for i = 1:length(numM)
    numM(i) = (m(i) - 5*log(distluz(omm, om, z(i))))/dm(i)^2;
end

denomM = zeros(1, 557);
for j=1:length(denomM)
    denomM(j) = 1/(dm(j))^2;
end

Mm = sum(numM)/sum(denomM);
disp(Mm);

% for j = 1:length(m)
%         dist(j) = distluz(0.3, -1, z(j));
% end

numchimarg = zeros(1, 557);
for i=1:length(numchimarg)
        numchimarg(i) = ((m(i) - 5*log(distluz(omm, om, z(i))) - Mm)^2)/(dm(i))^2;
end

chi = sum(numchimarg);
disp(chi);

