import nose,os
proteus=os.getenv('PROTEUS')
nose.run(defaultTest=os.path.join(proteus,'externalPackages','mpi4py','test'),argv=['runtests.py','-v','--exclude','test_spawn'])
print("PEXPECT_EXIT")
