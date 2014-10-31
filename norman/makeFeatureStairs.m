load trainsamples.mat

BIN_COUNT = 40;

for i=[1:size(Sample1.X,2)]
  bins1 = hist(Sample1.X(:,i), BIN_COUNT);
  bins2 = hist(Sample2.X(:,i), BIN_COUNT);
  bins3 = hist(Sample3.X(:,i), BIN_COUNT);

  figure();
  set(gcf,'Visible', 'off');

  stairs([1:BIN_COUNT], [(bins1/norm(bins1))', (bins2/norm(bins2))', (bins3/norm(bins3))'])
  legend('Sample1', 'Sample2', 'Sample3');
  title(sprintf('Feature %03d', i));

  saveas(gcf, sprintf('feature%03d.pdf', i));
end