import nose,os
proteus=os.getenv('PROTEUS')
nose.run(defaultTest=os.path.join(proteus,'externalPackages','mpi4py','test'))
print("PEXPECT_EXIT")
