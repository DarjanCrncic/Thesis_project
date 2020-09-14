
#include "main.h"
#include "functions.h"
#include <string.h>
#include <stdio.h>
#include "EventRecorder.h"              // Keil.ARM Compiler::Compiler:Event Recorder


#define dataSize 2000
#define crosPointNumber 1000

#define sizeOfDataAll 2974

#define wi 0.5
#define t 160
#define largeClustPerc 1 //0.8
#define mergePerc 0.8
#define maxdist 1024
#define maxClust 500
#define minMax 2
#define simMatrixSize 100

float data_freq_all[dataSize]={0};
float data_time_all[dataSize]={0};

float data_freq[dataSize]={0};
float data_time[dataSize]={0};

// array init.			
short int neighIndices[dataSize] = {0};
short int seeds[dataSize] = {0};
short int orderedList[dataSize] = {0};
float coreDistList[dataSize];
signed char procesList[dataSize] = {0};
float reachDistList[dataSize] = {0};
float orderedReachDistList[dataSize] = {0};

short int clusterIndices[dataSize] = {0};
short int prevVectorIndices[crosPointNumber]= {0};
short int currVectorIndices[crosPointNumber]= {0};
short int prevCluster[simMatrixSize]={0};
short int currCluster[simMatrixSize]={0};
short int newClustTest[simMatrixSize]={0};
	
int setClustersStart[maxClust]={0};
int setClustersEnd[maxClust]={0};
	
float timeMin[simMatrixSize]={0};
float timeMax[simMatrixSize]={0};
float yMinEnd[simMatrixSize]={0};
float yMaxEnd[simMatrixSize]={0};
float yMinStart[simMatrixSize]={0};
float yMaxStart[simMatrixSize]={0};
float prevyMaxEnd[simMatrixSize]={0};
float prevyMinEnd[simMatrixSize]={0};


UART_HandleTypeDef hlpuart1;

void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_LPUART1_UART_Init(void);

//retargeting for pritf and scanf
	struct __FILE
	{
		int dummyVar; //Just for the sake of redefining __FILE, not used
	};

	FILE __stdout;
	FILE __stdin;

	int fputc(int c, FILE * stream)
	{
	HAL_UART_Transmit(&hlpuart1,(uint8_t *)&c, 1, HAL_MAX_DELAY);
		return c; //return the character written to denote a successfull write
	}

	int fgetc(FILE * stream)
	{
		uint8_t c=0;
			//if (__HAL_UART_GET_FLAG(&hlpuart1, UART_FLAG_ORE))
			//	__HAL_UART_CLEAR_FLAG(&hlpuart1, UART_FLAG_ORE);
		HAL_UART_Receive(&hlpuart1,(uint8_t *)&c, 1, HAL_MAX_DELAY);
		return c;
	}
//////////////////////////////////////START/////////////////////////////////////////
// data initialization

int main(void)
{
	
  HAL_Init();
  SystemClock_Config();
  MX_GPIO_Init();
  MX_LPUART1_UART_Init();
	EventRecorderInitialize(EventRecordAll, 1U); //Keil Event Recorder used for timing measurment
	
	volatile int i, clustNum, newClustIndex=0, iter=0, reset, Nmin, eps, crossPoint, dataSizeTest, step=0;
	float timeStampStart, timeStampEnd;
	
	// data transport over UART
	
	HAL_UART_Init(&hlpuart1);	
	
	for(i=0;i<sizeOfDataAll;i++){
		scanf("%f",&data_freq_all[i]);
	}
	for(i=0;i<sizeOfDataAll;i++){
		scanf("%f",&data_time_all[i]);
	}
	
	while (1)
  {
	HAL_UART_Init(&hlpuart1);	
		
	// flag for software reset
	scanf("%d", &reset);
		
	scanf("%d", &dataSizeTest);	
	scanf("%d", &crossPoint);	
	scanf("%d", &eps);
	scanf("%d", &Nmin);
			
	// refreshing all arrays
	for(i=0;i<dataSizeTest;i++){
		
		data_freq[i]=0;
		data_time[i]=0;
		
		coreDistList[i]=-1;
		reachDistList[i]=-1;	
		neighIndices[i]=0;
		seeds[i]=0;
		clusterIndices[i]=0;
		procesList[i]=0;
		orderedList[i]=0;
	}	
	
	for(i=0;i<maxClust;i++){
		setClustersStart[i]=0;
		setClustersEnd[i]=0;
	}
	
	for(i=0;i<simMatrixSize;i++){
		newClustTest[i]=0;
		currCluster[i]=0;
		
		timeMin[i]=0;
		timeMax[i]=0;
		yMinEnd[i]=0;
		yMaxEnd[i]=0;
		yMinStart[i]=0;
		yMaxStart[i]=0;
	}
	
	// software reset
	if(reset==1){
		iter=0;
		newClustIndex=0;
		timeStampStart=0;
		step=0;
	}
	
	for(i=0;i<dataSizeTest;i++){
		data_freq[i]=data_freq_all[step*(dataSizeTest-crossPoint)+i];
		data_time[i]=data_time_all[step*(dataSizeTest-crossPoint)+i];
	}
	step=step+1;
 		
	// end initialization
	
	// optics algorithm, extracting reachability distance plot 
	
	//EventStartC(0); 	
	optics ( eps, Nmin, dataSizeTest );
	//EventStopC(0); 
	
	// extracting ordered reachability distance plot
	for(i=0;i<dataSizeTest;i++){
		orderedReachDistList[i]=reachDistList[orderedList[i]];
		
		if(orderedReachDistList[i]==-1)
			orderedReachDistList[i]=maxdist;
	}
	 		
	// gradient clustering algorithm
	//EventStartC(1); 	 	
	clustNum=gradientClustering( Nmin, t, wi, largeClustPerc, mergePerc, minMax, dataSizeTest);
	//EventStopC(1);
	
	// geting cluster indices
	getClusterIndices( clustNum, dataSizeTest );

	// window merging
	//EventStartC(2); 	 	
	newClustIndex = opticsMerging(  clustNum, iter, newClustIndex, dataSizeTest, crossPoint );
	//EventStopC(2);
	timeStampEnd=data_time[dataSizeTest-1];
	
	// feature compression
	//EventStartC(3); 		
	featureExtraction( clustNum, timeStampEnd, timeStampStart, Nmin, dataSizeTest );
  
	
	iter=1;	
	timeStampStart=timeStampEnd;
	
	for(i=0;i<simMatrixSize;i++){
		prevCluster[i]=currCluster[i];
		prevyMaxEnd[i]=yMaxEnd[i];
		prevyMinEnd[i]=yMinEnd[i];
	}
	for(i=0;i<crossPoint;i++){
		prevVectorIndices[i]=clusterIndices[dataSizeTest-crossPoint+i];
	}
	//EventStopC(3);
	
	//EventStopC(0);
	
	// sending results to MATLAB
	
	printf("%d\n\r", clustNum);
	
	for(i=0;i<dataSizeTest;i++){
		printf("%hd\n\r", clusterIndices[i]);
	}
	for(i=0;i<clustNum;i++){
		//printf("%d\n\r", setClustersStart[i]);
		//printf("%d\n\r", setClustersEnd[i]);	
		printf("%hd\n\r", currCluster[i]);
	}
	for(i=1;i<=clustNum;i++){
		printf("%f\n\r", timeMin[i]);
		printf("%f\n\r", timeMax[i]);
		printf("%f\n\r", yMinStart[i]);
		printf("%f\n\r", yMaxStart[i]);
		printf("%f\n\r", yMinEnd[i]);
		printf("%f\n\r", yMaxEnd[i]);
	}
//	for(i=0;i<dataSizeTest;i++){
//		printf("%f\n\r", orderedReachDistList[i]);
//	}

	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// clock, UART and PIN configuration
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};
  RCC_PeriphCLKInitTypeDef PeriphClkInit = {0};

  /** Configure the main internal regulator output voltage 
  */
  if (HAL_PWREx_ControlVoltageScaling(PWR_REGULATOR_VOLTAGE_SCALE1_BOOST) != HAL_OK)
  {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB busses clocks 
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 2;
  RCC_OscInitStruct.PLL.PLLN = 30;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = RCC_PLLQ_DIV2;
  RCC_OscInitStruct.PLL.PLLR = RCC_PLLR_DIV2;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB busses clocks 
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_5) != HAL_OK)
  {
    Error_Handler();
  }
  PeriphClkInit.PeriphClockSelection = RCC_PERIPHCLK_LPUART1;
  PeriphClkInit.Lpuart1ClockSelection = RCC_LPUART1CLKSOURCE_PCLK1;
  if (HAL_RCCEx_PeriphCLKConfig(&PeriphClkInit) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief LPUART1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_LPUART1_UART_Init(void)
{

  /* USER CODE BEGIN LPUART1_Init 0 */

  /* USER CODE END LPUART1_Init 0 */

  /* USER CODE BEGIN LPUART1_Init 1 */

  /* USER CODE END LPUART1_Init 1 */
  hlpuart1.Instance = LPUART1;
  hlpuart1.Init.BaudRate = 115200;
  hlpuart1.Init.WordLength = UART_WORDLENGTH_8B;
  hlpuart1.Init.StopBits = UART_STOPBITS_1;
  hlpuart1.Init.Parity = UART_PARITY_NONE;
  hlpuart1.Init.Mode = UART_MODE_TX_RX;
  hlpuart1.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  hlpuart1.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;
  hlpuart1.Init.ClockPrescaler = UART_PRESCALER_DIV1;
  hlpuart1.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;
  hlpuart1.FifoMode = UART_FIFOMODE_DISABLE;
  if (HAL_UART_Init(&hlpuart1) != HAL_OK)
  {
    Error_Handler();
  }
  if (HAL_UARTEx_SetTxFifoThreshold(&hlpuart1, UART_TXFIFO_THRESHOLD_1_8) != HAL_OK)
  {
    Error_Handler();
  }
  if (HAL_UARTEx_SetRxFifoThreshold(&hlpuart1, UART_RXFIFO_THRESHOLD_1_8) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN LPUART1_Init 2 */

  /* USER CODE END LPUART1_Init 2 */

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOG_CLK_ENABLE();
  HAL_PWREx_EnableVddIO2();

}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */

  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(char *file, uint32_t line)
{ 
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     tex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
