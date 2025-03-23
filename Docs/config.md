![](UNEDA%20stack.png)

The following files should not be primary targets in a makefile. Instead, they are included as part of the compile/link process. This is done for scoping reasons. Variables and code can have different scope. Apart from the traditional within block, within module and global, there is shared between a few modules via #include-ing the second (dependent) file from the first one.

In TCL:
+ TCLevalp.c
+ TCLwarp.c

In DTL/SML:
+ DTLautoscale.c
+ DTLdominance.c
+ SMLlayer.c

In CAR:
+ CARrank.c

For historical reasons, TCL is a separate library. It was developed before the others and was for a while the only existing platform. In modern times, all of those libraries can be linked together into one piece, forming a .dll (Win) or .so (Linux) for easy use by apps.
