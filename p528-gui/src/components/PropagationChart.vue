<template>
  <ScatterChartGenerator
    :chart-data="chartData"
    :chart-options="chartOptions"
  />
</template>

<script>
import { Scatter as ScatterChartGenerator } from 'vue-chartjs/legacy'
import {
  Chart as ChartJS,
  Title,
  Tooltip,
  Legend,
  LineElement,
  LinearScale,
  CategoryScale,
  PointElement
} from 'chart.js'

ChartJS.register(
  Title,
  Tooltip,
  Legend,
  LineElement,
  LinearScale,
  CategoryScale,
  PointElement
)
export default {
  name: 'LineChart',
  components: { ScatterChartGenerator },
  props: {
    chartData: {
        type: Object,
        required: true
      },
    chartOptions: {
      type: Object,
      required: true
    }
  },
  data() {
    return {
            
    };
  }
}
</script>
