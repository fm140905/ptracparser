# MCNP Binary PTRAC Parser
Timestamps saved in ASCII-formatted PTRAC file is significantly limited in precision, which is not usable for coincidence analysis. This code 
- parses binary PTRAC output file output by MCNP6.2 (`master` branch) and MCNPX v2.7.0 (`mcnpx2.7` branch)
- generates pulse trains with accurate time stamps based on particle histories.