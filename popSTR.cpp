int main_1(int argc, char * argv[]);
int main_2(int argc, char * argv[]);
int main_3(int argc, char * argv[]);
int main_4(int argc, char * argv[]);
int main_5(int argc, char * argv[]);

int main(int argc, char * argv[])
{
	string subcommand = argv[1];
    if (subcommand == "computeReadAttributes")
        return main_1(argc-1, argv+1);
    else
    {
    	if (subcommand == "computePnSlippage")
    		return main_2(argc-1, argv+1);
    	else
    	{
    		if (subcommand == "computePnSlippageDefault")
    			return main_3(argc-1, argv+1);
    		else
    		{
    			if (subcommand == "msGenotyper")
    				return main_4(argc-1, argv+1);
    			else
    			{
    				if (subcommand == "msGenotyperDefault")
    					return main_5(argc-1, argv+1);
    				else
    				{
    					std::cerr << "Unrecognized subcommand, options are: computeReadAttributes, computePnSlippage, computePnSlippageDefault, msGenotyper and msGenotyperDefault.\n";
    					return 1;
    				}
    			}
    		}
    	}
    }
}