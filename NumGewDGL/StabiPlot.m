%% Plot von Stabilitätsgebieten nach R(z)

[X1,X2]=meshgrid(linspace(-5,5));
Z=complex(X1,X2);
%P=1+(Z.^2)/2+Z; % A
P=1+Z./(1-Z./2+(Z.^2)./12); % B
%P=1+Z./(1-(Z./2)); % C


surf(X1,X2,sqrt(imag(P).^2+real(P).^2));hold on;
surf(X1,X2,ones(size(X1)),'FaceColor',[1 0 0]);
view(0,-90)
title('Stabilitätsgebiet S von R(z) in blau')
xlabel('Realteil von Z')
ylabel('Imaginärteil von Z')