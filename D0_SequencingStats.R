#reads per sample
data<-readRDS("./06-outputfiles/C_data.rds")
reads<-data.frame(sample.id=colnames(data)[2:ncol(data)],reads=colSums(data[,2:ncol(data)]))
reads
min1<-min(reads$reads)
min1
min2<-min(reads$reads[reads$reads!=min(reads$reads)])
min2
median<-median(reads$reads)
median

reads %>%
  ggplot(aes(x=reads))+
  geom_histogram(binwidth=5000)

reads %>%
  ggplot(aes(x=1,y=reads))+
  geom_jitter()+
  scale_y_log10()
reads %>%
  ggplot(aes(x=1,y=reads))+
  geom_boxplot()+
  scale_y_log10()

reads %>%
  arrange(reads) %>%
  ggplot(aes (x=1:nrow(.), y=reads))+
  geom_line()
reads %>%
  arrange(reads) %>%
  ggplot(aes (x=1:nrow(.), y=reads))+
  geom_line()+
  coord_cartesian(xlim=c(0,30),ylim=c(0,50000))
reads %>%
  arrange(reads) %>%
  ggplot(aes (x=1:nrow(.), y=reads))+
  geom_line()+
  coord_cartesian(xlim=c(0,10),ylim=c(0,20000))

saveRDS(reads, "./06-outputfiles/D0_reads.rds")

datalong<-pivot_longer(data, 2:ncol(data), names_to="sample", values_to="reads")
datalong %>%
  group_by(sample) %>%
  summarize(n_seqs = sum (reads),
            n_sings = sum(reads==1),
            goods = 100 * (1 - n_sings / n_seqs)) %>%
  ggplot(aes(x=n_seqs, y=goods))+
  geom_point()

datalong %>%
  group_by(sample) %>%
  summarize(n_seqs = sum (reads),
            n_sings = sum(reads>=1 & reads>=5),
            goods = 100 * (1 - n_sings / n_seqs)) %>%
  ggplot(aes(x=n_seqs, y=goods))+
  geom_point()

datalong %>%
  filter(reads>=1 & reads<=5)

datalong %>%
  ggplot(aes(x=1:nrow(.), y=reads))+
  geom_point()+
  scale_y_log10()
