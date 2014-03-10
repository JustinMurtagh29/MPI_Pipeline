a = [2 NaN	3 6;
7 NaN NaN 7;
3 NaN 2	3;
7 NaN 4	8;
3	3	3	5;
5	6	4	9;
8	7	3	8;
3	NaN 4	7;
4	7	4	8;
10	4	2	7;
4	NaN	7	6;
7	NaN	NaN	1;
5	NaN	7	8;
7	NaN	3	9;
4	NaN	8	10;
9	NaN	4	1;
6	NaN	7	9;
2	NaN	3	9;
6	NaN	7	8];
%%
x = 1:2:10;
hist(a,x);
xlabel('Vote (scale 1-10)');
ylabel('Frequency');
legend('B4B-frontal (19 votes total)', 'B4B-parietal (5 votes total)', 'B4B-temporal (17 votes total)', 'TRacer (17 votes total)');
set(gca, 'XTick', x);
set(gca, 'XTickLabel', {'1-2', '3-4' '5-6' '7-8' '9-10'});
%%
nanmean(a)
nanstd(a)