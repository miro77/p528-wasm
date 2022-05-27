<template>
  <v-app>
    <v-card flat color="transparent">
      <v-card-text>
        <v-row>
          <v-col class="pr-4">
            <v-slider
              label="Terminal 1 Height (m)"
              v-model="profile.h_1__user_units"
              class="align-center"
              max="80000"
              min="1.5"
              step="0.5"
              hide-details
            >
              <template v-slot:append>
                <v-text-field
                  v-model="profile.h_1__user_units"
                  class="mt-0 pt-0"
                  hide-details
                  single-line
                  type="number"
                  style="width: 60px"
                ></v-text-field>
              </template>
            </v-slider>
          </v-col>
        </v-row>
        <v-row>
          <v-col class="pr-4">
            <v-slider
              label="Terminal 2 Height (m)"
              v-model="profile.h_2__user_units"
              class="align-center"
              max="80000"
              min="1.5"
              step="0.5"
              hide-details
            >
              <template v-slot:append>
                <v-text-field
                  v-model="profile.h_2__user_units"
                  class="mt-0 pt-0"
                  hide-details
                  single-line
                  type="number"
                  style="width: 60px"
                ></v-text-field>
              </template>
            </v-slider>
          </v-col>
        </v-row>
        <v-row>
          <v-col class="pr-4">
            <v-slider
              label="Frequency (MHz)"
              v-model="profile.f__mhz"
              class="align-center"
              max="30000"
              min="100"
              step="1"
              hide-details
            >
              <template v-slot:append>
                <v-text-field
                  v-model="profile.f__mhz"
                  class="mt-0 pt-0"
                  hide-details
                  single-line
                  type="number"
                  style="width: 60px"
                ></v-text-field>
              </template>
            </v-slider>
          </v-col>
        </v-row>
        <v-row>
          <v-col class="pr-4">
            <v-slider
              label="Time Percentage (%)"
              v-model="profile.time"
              class="align-center"
              max="99"
              min="1"
              step="1"
              hide-details
            >
              <template v-slot:append>
                <v-text-field
                  v-model="profile.time"
                  class="mt-0 pt-0"
                  hide-details
                  single-line
                  type="number"
                  style="width: 60px"
                ></v-text-field>
              </template>
            </v-slider>
          </v-col>
        </v-row>
        <v-row>
          <v-col class="pr-4">
            <v-slider
              label="Points"
              v-model="profile.nbr_of_points"
              class="align-center"
              max="2000"
              min="10"
              step="1"
              hide-details
            >
              <template v-slot:append>
                <v-text-field
                  v-model="profile.nbr_of_points"
                  class="mt-0 pt-0"
                  hide-details
                  single-line
                  type="number"
                  style="width: 60px"
                ></v-text-field>
              </template>
            </v-slider>
          </v-col>
        </v-row>
      </v-card-text>
      <v-btn @click="callP528" :disabled="isCalculating">Calc P.528! </v-btn>
    </v-card>
    <propagation-chart
      :chartData="chartData"
      :chartOptions="options"
    ></propagation-chart>
  </v-app>
</template>

<script>
import PropagationChart from "./components/PropagationChart.vue";
import { calc } from "../worker-api";

export default {
  name: "App",

  components: {
    "propagation-chart": PropagationChart,
  },
  data: () => ({
    isCalculating: false,
    chartData: {
      datasets: [
        {
          label: "Data One",
          backgroundColor: "#f87979",
          data: [40, 39, 10, 40, 39, 80, 40],
          showLine: true,
          pointRadius: 0,
          borderColor: "rgba(200, 0, 0, 1)",
        },
      ],
    },
    datacollection: {
      datasets: [],
    },
    options: {
      scales: {
        y: {
          min: 0,
          max: 300,
          reverse: true,
          type: "linear",
          position: "bottom",
        },
      },
      legend: {
        display: false,
      },
      // parsing: {
      //   xAxisKey: 'x',
      //   yAxisKey: 'y'
      // },
      responsive: true,
      maintainAspectRatio: false,
    },
    profile: {
      h_1__user_units: 5,
      h_2__user_units: 10000,
      time: 50,
      Polarization: 0,
      f__mhz: 5000,
      nbr_of_points: 50,
    },
  }),
  methods: {
    async callP528() {
      this.isCalculating = true;
      // this.chartData.datasets[0].data = [];
      this.intputJson = JSON.stringify(this.profile);
      this.resultString = await calc(this.intputJson);
      var result = JSON.parse(this.resultString);
      const coords = result.d__km.map((el, index) => ({
        x: el,
        y: result.A__db[index],
      }));
      this.chartData.datasets[0].data = coords;
      this.isCalculating = false;
    },
  },
};
</script>
