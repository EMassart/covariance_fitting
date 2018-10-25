%% This code loads the files error_bezier_Oct6.mat, error_bezier_S.mat and Results_piecewise_geod.mat, and plots the results
% to generate figure 3 of the note associated to the Skype meeting of the
% 9 October

clear all; close all; clc;

load error_bezier_quotient_SD.mat;
error_SD_bezierQ = error_SD_rec;
sol_SD_bezierG = sol_SD_rec;

load error_bezier_section_SD.mat;
error_SD_bezierS = error_SD_rec;
sol_SD_bezierS = sol_SD_rec;

load Results_piecewise_geod.mat;
error_SD_piecewise = error_SD_rec;
sol_SD_piecewise = sol_SD_rec;

[train_data, params_train] = data_points_external_frame(r);      % training data
[test_data, params, patch] = data_points_all_but_external(r);      % training data
s = size(test_data);

error_mean_bezierG = mean(error_SD_bezierQ);
error_mean_bezierS = mean(error_SD_bezierS,2);
error_mean_piecewise = mean(error_SD_piecewise,2);

E = [error_SD_bezierQ; error_SD_bezierS; error_SD_piecewise];
E2 = E;
E3 = E;
d_id = zeros(1,s(1)*s(2));
d_mean = zeros(1,s(1)*s(2));


for i_test = 1:s(1)*s(2)
    
    % current test point
    C_ref = test_data{i_test}*test_data{i_test}';
    
    % normalization with respect to the identity
    d_id(i_test) = norm(C_ref - eye(3024),'fro')^2;
    E2(:,i_test) = E(:,i_test)./d_id(i_test);
    
    % the four corners of the associated patch (required for normalization 
    % with respect to the average square distance towards the orners of the patch)
    Y = zeros(3024,20,4);
    Y(:,:,1) = train_data{patch{i_test}(2),patch{i_test}(1)};
    Y(:,:,2) = train_data{patch{i_test}(2)+1,patch{i_test}(1)};
    Y(:,:,3) = train_data{patch{i_test}(2),patch{i_test}(1)+1};
    Y(:,:,4) = train_data{patch{i_test}(2)+1,patch{i_test}(1)+1};

    d = zeros(1,4);
    for j = 1:4
        d(j) = norm(C_ref - Y(:,:,j)*Y(:,:,j)','fro')^2;
    end
    d_mean(i_test) = mean(d);
    
    % normalized error with respect to the average square distance towards the orners of the patch
    E3(:,i_test) = E(:,i_test)./d_mean(i_test);
end

params_mat = cat(3,params{:});  % get the value of the parameters associated to the different test points.
theta = reshape(params_mat(1,1,:),1,s(1)*s(2));
magn = reshape(params_mat(1,2,:),1,s(1)*s(2));
params_train_mat = cat(3,params_train{:});  % get the value of the parameters associated to the different test points.
s_train = size(params_train_mat);
theta_train = reshape(params_train_mat(1,1,:),1,s_train(3));
magn_train = reshape(params_train_mat(1,2,:),1,s_train(3));

% plot the error for the Bézier with exact geodesics surface
figure; hold on;
stem3(theta(1:28), magn(1:28), 100*E3(1,1:28), 'Color', 'r');
stem3(theta(29:end), magn(29:end), 100*E3(1,29:end), 'Color', 'b');
stem3(theta_train, magn_train, zeros(size(theta_train)), 'xk');
ax = gca();
xlabel('$\theta$','Interpreter','Latex','Fontsize',14);
ylabel('$W$','Interpreter','Latex','Fontsize',14);
zlabel('$E_{\mathrm{N}}(\hat C(\theta_i, W_i))(\%)$','Interpreter','Latex','Fontsize',14);
ax.XTick = 0:0.5:4;
ax.YTick = 4:1.5:13;
grid on;

% 
% figure; hold on;
% stem3(theta(1:28), magn(1:28), 100*E3(5,1:28), 'Color', 'r');
% stem3(theta(29:end), magn(29:end), 100*E3(5,29:end), 'Color', 'b');
% stem3(theta_train, magn_train, zeros(size(theta_train)), 'xk');
% ax = gca();
% xlabel('$\theta$','Interpreter','Latex','Fontsize',14);
% ylabel('$W$','Interpreter','Latex','Fontsize',14);
% zlabel('$E_{\mathrm{N}2}(\hat C(\theta_i, W_i))(\%)$','Interpreter','Latex','Fontsize',14);
% ax.XTick = 0:0.5:4;
% ax.YTick = 4:1.5:13;
% grid on;




%% Uncomment the following lines to generate figure 4 of the note associated 
% to the Skype meeting of the 9 October

% i_test = 5;
% anchor1 = train_data{3,1};   %  W = 10, theta = 0
% anchor2 = train_data{3,2};      %  W = 10, theta = 1
% 
% t = linspace(0,1,100);
% orth = orth_pol(anchor1'*anchor2);
% C_ref = test_data{i_test}*test_data{i_test}';
% d = zeros(1,100);
% for i = 1:length(t)
%     Y = (1-t(i))*anchor1 + t(i)*anchor2*orth';
%     d(i) = norm(Y*Y' - C_ref, 'fro')^2;
% end
% 
% figure; 
% plot(t, d, '.-b');
% xlabel('$\theta$','Interpreter','Latex','Fontsize',14);
% ylabel('$||C(0.5,10) - \gamma_{C(0,10) \rightarrow C(1,10)}(t)||^2_{\mathrm{F}}$','Interpreter','Latex','Fontsize',14);
% 
%% Uncomment the following lines to generate figure 5 of the note associated 
% to the Skype meeting of the 9 October
% i_test = 7;
% anchor1 = train_data{4,1};   %  W = 4, theta = 0
% anchor2 = train_data{4,2};      %  W = 7, theta = 0
% 
% t = linspace(0,1,100);
% orth = orth_pol(anchor1'*anchor2);
% C_ref = test_data{i_test}*test_data{i_test}';
% d = zeros(1,100);
% for i = 1:length(t)
%     Y = (1-t(i))*anchor1 + t(i)*anchor2*orth';
%     d(i) = norm(Y*Y' - C_ref, 'fro')^2;
% end
% 
% figure;
% plot(t, d, '.-b');
% xlabel('$\theta$','Interpreter','Latex','Fontsize',14);
% ylabel('$||C(0.5,13) - \gamma_{C(0,13) \rightarrow C(1,13)}(t)||^2_{\mathrm{F}}$','Interpreter','Latex','Fontsize',14);
% axis([0 max(t) 0 max(d)+10])

%% Uncomment the following lines to generate figure 6 of the note associated 
% to the Skype meeting of the 9 October
% i_test = 29;
% anchor1 = train_data{1,1};   %  W = 4, theta = 0
% anchor2 = train_data{2,1};      %  W = 7, theta = 0
% 
% t = linspace(0,1,100);
% orth = orth_pol(anchor1'*anchor2);
% C_ref = test_data{i_test}*test_data{i_test}';
% d = zeros(1,100);
% for i = 1:length(t)
%     Y = (1-t(i))*anchor1 + t(i)*anchor2*orth';
%     d(i) = norm(Y*Y' - C_ref, 'fro')^2;
% end
% 
% figure; 
% plot(t, d, '.-b');
% xlabel('$W$','Interpreter','Latex','Fontsize',14);
% ylabel('$||C(0,5.5) - \gamma_{C(0,4) \rightarrow C(0,7)}(t)||^2_{\mathrm{F}}$','Interpreter','Latex','Fontsize',14);
% 

%% Uncomment the following lines to generate figure 7 of the note associated 
% to the Skype meeting of the 9 October
% i_test = 2;
% anchor1 = test_data{1};   %  W = 4, theta = 0
% anchor2 = test_data{3};      %  W = 7, theta = 0
% 
% t = linspace(0,1,100);
% orth = orth_pol(anchor1'*anchor2);
% C_ref = test_data{i_test}*test_data{i_test}';
% d = zeros(1,100);
% for i = 1:length(t)
%     Y = (1-t(i))*anchor1 + t(i)*anchor2*orth';
%     d(i) = norm(Y*Y' - C_ref, 'fro')^2;
% end
% 
% figure;
% plot(t, d, '.-b');
% xlabel('$W$','Interpreter','Latex','Fontsize',14);
% ylabel('$||C(0.5,5.5) - \gamma_{C(0.5,4) \rightarrow C(0.5,7)}(t)||^2_{\mathrm{F}}$','Interpreter','Latex','Fontsize',14);

