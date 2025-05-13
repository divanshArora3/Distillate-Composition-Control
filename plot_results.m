time  = data.Time;
results = data.Data;
distillate_sp = results(:,1);
distillate = results(:,2);
reflux_ratio = results(:,3);
feed = results(:,4);
x_feed = results(:,5);

figure(1)
hold off

subplot(3,1,1)
hold off
plot(time,reflux_ratio,'b-','LineWidth',2)
axis([min(time) max(time) 0 11]);
legend('Reflux Ratio')
ylabel('Reflux Ratio')

subplot(3,1,2)
hold off
plot(time,x_feed,'g--','LineWidth',2)
hold on
plot(time,feed,'k-','LineWidth',2)
legend('Feed Composition','Feed Flow')
axis([min(time) max(time) 0.0 1.0]);
ylabel('Disturbances')

subplot(3,1,3)
hold off
plot(time,distillate_sp,'m-','LineWidth',2)
hold on
plot(time,distillate,'b:','LineWidth',2)
legend('Distillate SP','Distillate')
axis([min(time) max(time) ...
    min(min(distillate),min(distillate_sp))-0.1 ...
    1.0]);
ylabel('Distillate Conc')
xlabel('Time (min)')

% save data to text file
export = [time results];

save -ascii 'data.txt' export