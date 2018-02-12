library(ggplot2)
ggplot(cand_structs2, aes(Resolution)) + geom_density(adjust = 2/3, fill="skyblue", color="skyblue", alpha=0.3) +
  labs(x="Resolution in Angstroms",
       y="Count of Structures",
       title="Resolution of Filtered RCSB Structures w mRNA, S12")
ggsave(filename="./figures/resolution.png")
ggplot(cand_structs2, aes(mRNA_Length)) + geom_density(adjust = 1/3, fill="red", color="red", alpha=0.3) +
  labs(x="Length of Total mRNA (nt)",
       y="Count",
       title="Amount of mRNA in Filtered RCSB Structures w mRNA, S12")
ggsave(filename="./figures/mRNALength.png")

