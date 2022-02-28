# MCNP Binary PTRAC Parser
Timestamps saved in ASCII-formatted PTRAC file is significantly limited in precision, which is not usable for coincidence analysis. This code 
- parses binary PTRAC output file output by MCNPX2.7.0
- generates pulse trains with accurate time stamps based on particle histories.