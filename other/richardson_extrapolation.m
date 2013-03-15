% results of 1d CV simulation
% first column is h, second is dt, third is CV
a = [0.05, 0.01, 0.0736375371
0.025, 0.005, 0.0798080538
0.0125, 0.0025, 0.0770119017
0.00625, 0.00125, 0.0752304019
0.003125, 0.000625, 0.0746756238
0.0015625, 0.0003125, 0.0745052401
0.00078125, 0.00015625, 0.0744549664
0.000390625, 7.8125e-05, 0.074436781
0.0001953125, 3.90625e-05, 0.0744289901];

%% results of 1d monodomain model problem
%% first column is h, second the computed ||V(1)||_{Linf(L2)}
%% If uncomment this, change the for i=4:size(a,1) below to i=3:size(a,1)
%% and change CVh to a(1:end,2)
%a = [0.1, 0.95326;
%0.05, 0.987993;
%0.025, 0.996975;
%0.0125, 0.999242;
%0.00625, 0.999811;
%0.003125, 0.999953];

h = a(1:end,1);
dt = a(1:end,2);
CVh = a(1:end,3);

p = 0*h;
CVr = 0*CVh;
GCI = 0*CVh;

for i=4:size(a,1)
  w1 = CVh(i);
  w2 = CVh(i-1);
  w3 = CVh(i-2);

  h2_over_h1 = h(i-1)/h(i);

  pp = log( (w3-w2)/(w2-w1) ) / log(h2_over_h1);

  p(i) = pp;
  CVr(i) = w1 + (w1-w2)/(h2_over_h1^pp - 1);
  GCI(i) = abs((w1-w2)/w1) * 1.25/ (h2_over_h1^pp - 1);
end;


disp('h dt CVh CVr p GCI');
ret = [h dt CVh CVr p GCI]


