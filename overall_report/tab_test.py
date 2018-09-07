from tabulate import tabulate

columns=['Interface', 'IP', 'Status', 'Protocol']
sh_ip_int_br = [('FastEthernet0/0', '15.0.15.1', 'up', 'up'),
                ('FastEthernet0/1', '10.0.12.1', 'up', 'up'),
                ('FastEthernet0/2', '10.0.13.1', 'up', 'up'),
                ('Loopback0', '10.1.1.1', 'up', 'up'),
                ('Loopback100', '100.0.0.1', 'up', 'up')]

table = tabulate(sh_ip_int_br, headers=columns, tablefmt="latex", floatfmt=".2f")
print(table)
