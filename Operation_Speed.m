% Chun Yii Liu, Fast Algorithm, hw1
j = 1; k = 2; a = 1.234; b=2.345; N=1000000000;

% For empty loop --------------------
tic;
for i=1:N
    
end
T0 = toc;
fprintf('Empty loop:%f\n',T0);

% For + operator ---------------------
tic;
for i=1:N
    j = j+1;
    %a = a+1;
    %a = a+1.234;
end
T1 = toc;
fprintf('+ time + loop:%f\n',T1);
fprintf('Real + time:%f\n\n',T1-T0);

% For * operator ------------------------
tic;
for i=1:N
    k = k*2;
    %b = b*2;
    %b = b*1.234;
end
T1 = toc;
fprintf('* time + loop:%f\n',T1);
fprintf('Real * time:%f\n\n',T1-T0);

% For sin operator ------------------------
tic;
for i=1:N
    a = sin(a);
end
T1 = toc;
fprintf('sin time + loop:%f\n',T1);
fprintf('Real (sin) time:%f\n\n', T1-T0);

return