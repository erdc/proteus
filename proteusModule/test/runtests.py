import nose,os
proteus=os.getenv('PROTEUS')
nose.run(defaultTest=os.path.join(proteus,'proteusModule','test'),argv=['runtests.py','-v'])
print("PEXPECT_EXIT")
