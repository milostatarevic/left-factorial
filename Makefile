all:
	@$(MAKE) -C single
	@$(MAKE) -C range

clean:
	@$(MAKE) -C single clean
	@$(MAKE) -C range clean
