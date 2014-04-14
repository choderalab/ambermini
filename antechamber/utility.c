int spaceline(char *str) {
/*spaceline return 1 if str has only space or control characters*/
int i; 
int flag = 1;
	for(i=0;i<strlen(str);i++)
		if(isgraph(str[i])) {
			flag = 0;
			break;
		}
	return flag;
}
void newitoa(long num, char *str)
{
	int i = 0;
	int j;
	int type = 0;
	long tmpint, mod;
	char tmpchar[MAXCHAR];
	if(num < 0) {
		num *= -1;	
		type = 1;
	}
	strcpy(str, "");
	while (num >= 10) {
		tmpint = num / 10;
		mod = num - tmpint * 10;
		switch (mod) {
		case 0:
			tmpchar[i] = '0';
			break;
		case 1:
			tmpchar[i] = '1';
			break;
		case 2:
			tmpchar[i] = '2';
			break;
		case 3:
			tmpchar[i] = '3';
			break;
		case 4:
			tmpchar[i] = '4';
			break;
		case 5:
			tmpchar[i] = '5';
			break;
		case 6:
			tmpchar[i] = '6';
			break;
		case 7:
			tmpchar[i] = '7';
			break;
		case 8:
			tmpchar[i] = '8';
			break;
		case 9:
			tmpchar[i] = '9';
			break;
		}
		i++;
		num = tmpint;
	}
	switch (num) {
	case 0:
		tmpchar[i] = '0';
		break;
	case 1:
		tmpchar[i] = '1';
		break;
	case 2:
		tmpchar[i] = '2';
		break;
	case 3:
		tmpchar[i] = '3';
		break;
	case 4:
		tmpchar[i] = '4';
		break;
	case 5:
		tmpchar[i] = '5';
		break;
	case 6:
		tmpchar[i] = '6';
		break;
	case 7:
		tmpchar[i] = '7';
		break;
	case 8:
		tmpchar[i] = '8';
		break;
	case 9:
		tmpchar[i] = '9';
		break;
	}
	i++;
	if(type == 0) {
		for (j = 0; j < i; j++)
			str[j] = tmpchar[i - j - 1];
		str[i] = '\0';
	}
	else {
		str[0]='-';
		for (j = 0; j < i; j++)
			str[j+1] = tmpchar[i - j - 1];
		str[i+1] = '\0';
	}
	/*   printf("\n %5d",i);
	   printf("\n %s", tmpchar);
	   printf("\n %s", str); */
}

/*
void  dtoa (double num, char *str, int digital_num) {
char tmpchar1[MAXCHAR];
char tmpchar2[MAXCHAR];
char tmpchar3[MAXCHAR];
double a,b;
int c;
int count; 
	if(num==0) {
		strcpy(str, "0.");
		printf("\n%s", str); 
		return;
	}
	if(num<0) {
		num*=-1;
		strcpy(str, "-");
	}
	else
		strcpy(str,"");
	a = floor(num) ;	
	newitoa(a, tmpchar1);
	b = num -a ;
	strcpy(tmpchar2, "");
	tmpchar2[0]='\0';
	count = 0;
	while(b<1) {
		b*=10;
		if(count != 0) tmpchar2[count-1] = '0';
		count++;
	}
	tmpchar2[count-1] = '\0';
	count = 0;
	while(count < digital_num) {
		c = floor(b);
		switch(c) {
        		case 0:
                		tmpchar3[count] = '0';
                	break;
        		case 1:
                		tmpchar3[count] = '1';
                	break;
        		case 2:
                		tmpchar3[count] = '2';
                	break;
        		case 3:
                		tmpchar3[count] = '3';
                	break;
        		case 4:
                		tmpchar3[count] = '4';
                	break;
        		case 5:
                		tmpchar3[count] = '5';
                	break;
        		case 6:
                		tmpchar3[count] = '6';
                	break;
        		case 7:
                		tmpchar3[count] = '7';
                	break;
        		case 8:
                		tmpchar3[count] = '8';
                	break;
        		case 9:
                		tmpchar3[count] = '9';
                	break;
        	}
		count++;
		b=(b-c)*10;
	}
	tmpchar3[count] = '\0';
	strcat(str, tmpchar1);
	strcat(str,".");
	strcat(str, tmpchar2);
	strcat(str, tmpchar3);
}
*/
