library(ggplot2)
library(dplyr)

# Plot the resolution of all candidate structures
ggplot(cand_structs2, aes(Resolution)) + geom_density(adjust = 2/3, fill="skyblue", color="skyblue", alpha=0.3) +
  labs(x="Resolution in Angstroms",
       y="Count of Structures",
       title="Resolution of Filtered RCSB Structures w mRNA, S12")
ggsave(filename="./figures/resolution.png")

# Plot how much mRNA is available across structures
ggplot(cand_structs2, aes(mRNA_Length)) + geom_density(adjust = 1/3, fill="red", color="red", alpha=0.3) +
  labs(x="Length of Total mRNA (nt)",
       y="Count",
       title="Amount of mRNA in Filtered RCSB Structures w mRNA, S12")
ggsave(filename="./figures/mRNALength.png")

MethodData <- cand_structs2 %>%
  group_by(Method) %>%
  summarize(count = n()) %>%
  mutate(pct = count / sum(count))

# Plot how the various candidate structures are imaged
ggplot(MethodData, aes(x=Method, fill=Method, y=pct)) + geom_bar(stat="identity") + 
  labs(y="Count of Structures", 
       title="How the Candidate Filtered RCSB Structures are Imaged") +
  geom_text(aes(label=paste("n=", scales::comma(count), sep=""), vjust=2), color="white")
ggsave(filename="./figures/imaging.png")

library(stringr)
# Want to see the nucleotide Distributions for the Structures
A_count <- sum(str_count(cand_structs2w$mRNA_Sequences, "A"))
U_count <- sum(str_count(cand_structs2w$mRNA_Sequences, "U"))
C_count <- sum(str_count(cand_structs2w$mRNA_Sequences, "C"))
G_count <- sum(str_count(cand_structs2w$mRNA_Sequences, "G"))

counts = c(A_count, U_count, C_count, G_count)
total_count <- sum(counts)
freq_counts <- sapply(counts, function(count) count/total_count)
labels <- c("A","U","C","G")

df <- data.frame(labels, freq_counts, counts)

ggplot(df, aes(x=labels, y=freq_counts, fill=labels)) + geom_bar(stat="identity") + 
  labs(y="Frequency of Nucleotide", 
       title="Nucleotide Frequency in Structures") +
  geom_text(aes(label=paste("n=", scales::comma(counts), sep=""), vjust=2), color="white")

ggsave("./figures/nucleotide_frequencies.png")

# df <- df%>%
#   dplyr::arrange(counts) %>%
#   mutate(labels=factor(labels, levels =labels))
# 
# ggplot(df, aes(x=labels, y=freq_counts, fill=labels)) + geom_bar(stat="identity") + 
#   labs(y="Frequency of Nucleotide", 
#        title="Nucleotide Frequency in Structures") +
#   geom_text(aes(label=paste("n=", scales::comma(counts), sep=""), vjust=2), color="white")
# 
# ggsave("./figures/nucleotide_frequencies_ordered.png")
