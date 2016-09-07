//----------------------------------------------------------------------------
// Utility functions
//----------------------------------------------------------------------------
function ToRadians(degrees: number):number {
    return (Math.PI * degrees) / 180.0;
}

function ToDegrees(rads:number):number{
    return rads * 180.0/ Math.PI;
}

function Rand(min: number, max: number): number {
    return Math.random() * (max - min) + min;
}

function RandNoise(mu:number,sigma:number):number{
    return Math.random() * sigma + mu;
}

//ref https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/random
// Returns a random integer between min (included) and max (included)
// Using Math.round() will give you a non-uniform distribution!
function RandomIntInclusive(min, max):number {
  return Math.floor(Math.random() * (max - min + 1)) + min;
}

//ref https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/random
// Returns a random integer between min (included) and max (excluded)
// Using Math.round() will give you a non-uniform distribution!
function RandomInt(min, max):number {
  return Math.floor(Math.random() * (max - min)) + min;
}

function Clamp(min: number, max: number, val: number):number {
    return (val <= min ? min : (val >= max ? max : val));
}

function leftPad(v:number,base:number):string{
    return v<10? ('0'+v.toString(base)) : (v.toString(base));
}

function RGB2Hex(r:number,g:number,b:number):string{
    var rs = leftPad(Clamp(0,255,r), 16);
    var gs = leftPad(Clamp(0,255,g), 16);
    var bs = leftPad(Clamp(0,255,b), 16);
    return '#'+rs+gs+bs;
}

function ClampedSub(v:number,k:number):number{
    v -= (v-k >=0) ? k:0
    return v;
}

//exponential scaling to a range
//For example scale v to exp range (0.1,10), which means scale
// v to a linear scale of log(0.1)=-1,log(10)=1 and then taking exponent of the resulting value
function ExpScale(min:number,max:number,v:number):number{
    var min_e = Math.log(min);
    var max_e = Math.log(max);
    return Math.exp(v*(max_e-min_e) + min_e)
}

//scale a value of range[0,1] linearly to the range (min,max)
function LinScale(min:number,max:number,v:number):number{
    return v * (max-min) + min;
}

//This function returns an array of indices randomly shuffled.
function ShuffledIndices(length:number):Array<number>{
    var indices = new Array<number>(length);
    var retval = new Array<number>(length);
    for(var i=0;i<length;i++){
        indices[i] = 0;
    }

    var j = 0;
    while(j < length){
        var i = RandomInt(0,length);
        if(indices[i]==0){
            indices[i] = 1;
            retval[j++] = i;
        }
    }
    return  retval;
}

function GetRandomFloatArray(length:number,min:number,max:number) {
    var retval = new Array<number>(length);
    for(var i=0;i<length;i++){
        retval[i] = Rand(min,max);
    }
    return retval;
}

function GetRandomIntArray(length:number,min:number,max:number) {
    var retval = new Array<number>(length);
    for(var i=0;i<length;i++){
        retval[i] = RandomInt(min,max);
    }
    return retval;
}

//Generates a random cartesian coordinate inside a given polar components
function RandomCartesianCoordinate(r:number,theta:number=Math.PI*2):Vector2D{
    var t = Rand(0,theta);
    var x = r * Math.cos(t);
    var y = r * Math.sin(t);
    return new Vector2D(x,y);
}
//----------------------------------------------------------------------------
// Vector classes
//----------------------------------------------------------------------------
class Vector2D{
    constructor(public x:number,public y:number){
    }
    public norm(){
        return Math.sqrt( this.x *this.x + this.y*this.y);
    }
    public add(v:Vector2D):Vector2D{
        return new Vector2D(this.x + v.x, this.y + v.y);
    }
    public sub(v:Vector2D):Vector2D{
        return new Vector2D(this.x - v.x, this.y - v.y);
    }
    public dot(v:Vector2D):number{
        return this.x * v.x + this.y * v.y;
    }
    public angle(v:Vector2D):number{
        var j = this.dot(v);
        var k = (this.norm() * v.norm());
        var d = j / k;
        // console.log('angle',d,'dot',this.dot(v),'|V||V|',k);
        return Math.acos(d);
    }
}

//----------------------------------------------------------------------------
// Simple Chart class
//----------------------------------------------------------------------------
enum ChartType{Linear,Log}
class Chart{
    private Values:Array<number> = new Array<number>();
    private margin = 5;
    private min:number = Number.MAX_VALUE;
    private max:number = Number.MIN_VALUE;
    private v_range:number = 1;
    private h_range:number = 1;
    private y_r:number = 1;
    private x_r:number = 1;
    constructor(public ctx:CanvasRenderingContext2D, public position:Vector2D,
                public size:Vector2D,public Type:ChartType){
    }

    public Clear(){
        this.Values = new Array<number>();
        this.min = Number.MAX_VALUE;
        this.max = Number.MIN_VALUE;
        this.v_range = 1;
        this.h_range = 1;
        this.y_r = 1;
        this.x_r = 1;
    }

    public Update(value:number){
        if(this.Type == ChartType.Log)
            if(value ==0)
                value = 0.000000001; //THIS IS A HACK TO PREVENT -inifinity values
        this.Values.push(value);

        if(value<=this.min)this.min=value;
        if(value>=this.max)this.max=value;

        this.v_range = this.max-this.min;
        if(this.Type == ChartType.Log)
            this.v_range = Math.log(this.max)-Math.log(this.min);

        this.h_range = this.Values.length;
        this.y_r = (this.size.y - (this.margin*2)) / this.v_range;
        this.x_r = (this.size.x - (this.margin*2)) / this.h_range;

        this.Draw();
    }

    public Draw(){
        if(this.Values.length==0)
            return;

        if(this.min==this.max)
            return;
        this.ctx.save();
        this.ctx.setTransform(1, 0, 0, 1, 0, 0);
        this.ctx.clearRect(this.position.x,this.position.y,this.size.x,this.size.y);
        this.ctx.lineWidth = 4;
        this.ctx.strokeStyle = 'green'
        this.ctx.strokeRect(this.position.x,this.position.y,this.size.x,this.size.y);
        // this.ctx.fillStyle = 'gainsboro';
        // this.ctx.fillRect(this.position.x,this.position.y,this.size.x,this.size.y);

        var v_offset = this.position.y + this.size.y
        var x = this.position.x + this.margin;

        var y = v_offset - ((this.Values[0]-this.min) * this.y_r + this.margin);
        if(this.Type == ChartType.Log)
            y = v_offset - ((Math.log(this.Values[0])-Math.log(this.min)) * this.y_r + this.margin);

        this.ctx.beginPath();
        this.ctx.moveTo(x,y);
        this.ctx.lineWidth = 1;
        this.ctx.strokeStyle = 'red';
        var xx = 0;
        if(this.Type == ChartType.Log){
            for(var v of this.Values){
                x = this.position.x + xx * this.x_r + this.margin
                y = v_offset - ((Math.log(v) - Math.log(this.min)) * this.y_r + this.margin)
                this.ctx.lineTo(x,y);
                this.ctx.stroke();
                xx++;
            }
        }else{
            for(var v of this.Values){
                x = this.position.x + xx * this.x_r + this.margin
                y = v_offset - ((v - this.min) * this.y_r + this.margin)
                this.ctx.lineTo(x,y);
                this.ctx.stroke();
                xx++;
            }
        }

        this.ctx.closePath();
        //draw axis markers

        //end marker
        var xx =(this.Values.length-1);
        this.draw_horiz_text(''+xx,this.size.x - this.margin);
        //mid marker
        xx = Math.round(this.Values.length/2);
        this.draw_horiz_text(''+xx,this.size.x/2 - this.margin);
        this.ctx.restore();

        //max val markert
        var v= this.max;
        var yy = this.margin*4;
        this.draw_vert_text(''+v.toFixed(2),yy);

        var v= this.min;
        var yy = this.size.y + this.margin;
        this.draw_vert_text(''+v.toFixed(2),yy);
    }

    private draw_horiz_text(s:string,xx:number){
        var x = this.position.x + xx;
        var y = this.position.y + this.size.y + this.margin
        this.ctx.font = "40px arial";
        this.ctx.fillStyle='red'
        this.ctx.fillText(s,x,y);
        this.ctx.fillRect(x,y,3,10);
    }

    private draw_vert_text(s:string,yy:number){
        var x = this.position.x - this.margin;
        var y = this.position.y + yy
        this.ctx.font = "40px arial";
        this.ctx.fillStyle='red'
        this.ctx.fillText(s,x,y);
        this.ctx.fillRect(x,y,10,3);
    }

}

//----------------------------------------------------------------------------
// CTRNN class
//----------------------------------------------------------------------------
class Network{
    public Z:Array<number>;
    public Y:Array<number>;
    public I:Array<number>;

    private RuleFunctors:{(i:number,j:number):number;}[]

    constructor(
                public W:Array<Array<number>>,
                public B:Array<number>,
                public Tau:Array<number>,
                public PlasticRules:Array<Array<number>>,
                public Eta:Array<Array<number>>
                ){
                    this.Z = new Array<number>(B.length);
                    this.Y = new Array<number>(B.length);
                    this.I = new Array<number>(B.length);
                    for(var j=0;j<this.Z.length;j++){
                        this.Z[j] = 0;
                        this.Y[j] = 0;
                        this.I[j] = 0;
                    }

                    // this.RuleFunctors = [this.Rule0,this.Rule1,this.Rule2,this.Rule3];
                }
    private ComputeActivations(){
        for(var j=0;j<this.Y.length;j++){
            this.Z[j] = 1 / (1 + Math.exp(-this.Y[j]-this.B[j]));
        }
        return this.Z;
    }

    private ComputePlasticFacilitation(j:number){
        var yj = this.Y[j]
        var bj = this.B[j];
        var t = yj+bj;
        if (t < -4)
            return -1;
        else if(t>4)
            return +1;
        else if(t > -2 && t < 2)
            return 0;
        else if(t > -4 && t < -2)
            return 0.5 * (t+4) - 1;
        else if(t > +2 && t < +4)
            return 0.5 * (t-2);
    }
    // private Rule0(i:number,j:number):number{
    //     return this.Eta[i][j] * this.ComputePlasticFacilitation(j) * this.Z[i] * this.Z[j];
    // }
    // private Rule1(i:number,j:number):number{
    //     var Zo = (this.W[i][j] + 8) / 16;
    //     return this.Eta[i][j] * this.ComputePlasticFacilitation(j) * (this.Z[i] - Zo) * this.Z[j];
    // }
    // private Rule2(i:number,j:number):number{
    //     var Zo = (this.W[i][j] + 8) / 16;
    //     return this.Eta[i][j] * this.ComputePlasticFacilitation(j) * this.Z[i] * (this.Z[j] - Zo);
    // }
    // private Rule3(i:number,j:number):number{
    //     return 0;
    // }

    private ApplyPlasticRule(rule:number,i:number,j:number):number{

        switch (rule) {
            case 0: {
                return this.Eta[i][j] * this.ComputePlasticFacilitation(j) * this.Z[i] * this.Z[j];
            }
            case 1: {
                var Zo = (this.W[i][j] + 8) / 16;
                return this.Eta[i][j] * this.ComputePlasticFacilitation(j) * (this.Z[i] - Zo) * this.Z[j];
            }
            case 2: {
                var Zo = (this.W[i][j] + 8) / 16;
                return this.Eta[i][j] * this.ComputePlasticFacilitation(j) * this.Z[i] * (this.Z[j] - Zo);
            }
            case 3: return 0;
            default:
                break;
        }
    }

    private ApplyPlasticRules():number{
        var homeostatic_neuronal_ratio = 0; //number of neurons that behave without plasticity
        for(var i=0;i<this.Y.length;i++){
            for(var j=0;j<this.Y.length;j++){
                var rule = this.PlasticRules[i][j];
                //var Wd = this.RuleFunctors[rule](i,j);
                var Wd = this.ApplyPlasticRule(rule,i,j);
                var delta = (8-Math.abs(this.W[i][j]));
                var new_val = this.W[i][j] + Wd * delta
                this.W[i][j] = Clamp(-8,8,new_val)

                //homeostatic_neuronal_ratio += ((Wd > 0 || Wd < 0) ? 0:1);
                homeostatic_neuronal_ratio += ((Math.abs(Wd) < 1e-10) ? 1:0);
            }
        }
        return homeostatic_neuronal_ratio / (this.B.length * this.B.length); //normalized to number of connections
    }

    public Evaluate():number{
        var homeostatic_neuronal_ratio = this.ApplyPlasticRules();

        for(var i=0;i<this.Y.length;i++){
            var ev = 0;
            for(var j=0;j<this.Y.length;j++){
                ev += (this.W[i][j] * this.Z[j]);
            }
            var dy_dt = 1/this.Tau[i] * (-this.Y[i] + ev + this.I[i]);

            //Forward Euler Integration
            var h = 0.1;
            this.Y[i] = this.Y[i] + h * dy_dt;
        }
        this.ComputeActivations();
        return homeostatic_neuronal_ratio;
    }
}

//----------------------------------------------------------------------------
// Genetic Algorithm class
//----------------------------------------------------------------------------
class SteadyStateGeneticAlgorithm{
    public Population:Array<Array<number>>;
    public Fitness:Array<number>;
    public CurrentGeneration:number = 0;

    constructor(public PopulationCount:number, public pMutation:number,
                                                public pCrossover:number, private FloatChromosomeLength:number,private IntChromosomeLength:number){
        this.Population = new Array<Array<number>>(PopulationCount);
        this.Fitness = new Array<number>(PopulationCount);
        var i=0;
        //Initialize
        //NOTE: Chromosome has 64 element int array appended for plasticity rules
        while(i < PopulationCount){
            this.Population[i] = GetRandomFloatArray(FloatChromosomeLength,0,1).concat(GetRandomIntArray(this.IntChromosomeLength,0,4));
            this.Fitness[i] = 0;
            i++;
        }
    }

    private Mutate(winner:number,loser:number,pm:number,pc:number){
        var len = this.FloatChromosomeLength + 64;
        var p_mutate = GetRandomFloatArray(len,0,1);
        var mutations = GetRandomFloatArray(len-64,-0.01,0.01).concat(GetRandomIntArray(64,-1,2));
        var p_crossover = GetRandomFloatArray(len,0,1);

        for(var i=0;i<len;i++){
            //crossover
            if(p_crossover[i]>pc){
                this.Population[loser][i] = this.Population[winner][i];
            }
            //mutate
            if(p_mutate[i]>pm){
                var val = this.Population[loser][i] + mutations[i];

                if(i<this.FloatChromosomeLength){
                    if(val>=0 && val <=1)  //float values in the range [0,1]
                        this.Population[loser][i] = val
                }else{
                    if(val>=0 && val <=3) //int values in the range [0,4]
                        this.Population[loser][i] = val
                }
            }
        }
    }

    public Evolve(fitness:Array<number>){
        this.Fitness = fitness; //TODO: is this necessary, may be just use as is and remove class member

        // steady state GA
        var shuffled_indices = ShuffledIndices(this.PopulationCount);
        //Tournament selection between two random individuals with no duplicates
        for(var k=0;k<this.PopulationCount/2;k++){
            var I1= shuffled_indices[2*k];
            var I2= shuffled_indices[2*k+1];

            var winner = I1, loser = I2; // assume these values
            if(this.Fitness[I1] <= this.Fitness[I2]){//maximise fitness
                winner = I2;
                loser = I1;
            }
            //  mututate and apply crossover to the losers
            this.Mutate(winner,loser,this.pMutation,this.pCrossover);
            // console.log('winner,loser',winner,loser)
        }
        this.CurrentGeneration++;
    }

    public SerializeTo():string{
        return (JSON.stringify(this.Population)+','+JSON.stringify(this.Fitness));
    }
    public Deserialize(s:string){
        var i1  = s.indexOf(']],')
        var pop = s.substr(0,i1+2);
        var fit = s.substr(i1+3)
        this.Population = JSON.parse(pop);
        this.Fitness = JSON.parse(fit);
    }
}

//============================================================================
// Simulation definition classes
//============================================================================
interface ISituatedObject {
    Id: number; //Unique identifier for this object
    Position:Vector2D;
    Radius: number;
    //update internal state
    //width & height = canvas width & height
    Update(vect:Vector2D,r:number);
    Draw(ctx: CanvasRenderingContext2D);
}
//----------------------------------------------------------------------------
// Bot class
//----------------------------------------------------------------------------
class Bot implements ISituatedObject{
    public Color = 'green';
    public Theta:number = 0.0;
    public SensorAngle:number = 45.0;
    public LeftMotorSpeed:number = 0.0;
    public RightMotorSpeed:number = 0.0;
    public LeftMotorGain:number = 0.0;
    public RightMotorGain:number = 0.0;
    public LeftSensorGain:number = 0.0;
    public RightSensorGain:number = 0.0;

    public IsInverted:boolean = false; //denotes if the sensors are inverted
    //FITNESS METRICS-------
    public SolarVicinityTime:number = 0.0;
    public IntegralSolarDistanceOverTime = 0.0;
    public HomeostaticNeuronalRatio = 0.0;
    public CurrentPlasticNeuronalRatio = 0.0;
    public Di = 0.0; //initial distance to the source
    public Df = 0.0; //final distance to the source
    //------------
    public Ctrnn:Network;

    constructor(public Id:number, public Position:Vector2D,
        public Radius:number, angle:number, ctrnn:Network){
            this.Theta = angle;
            this.Ctrnn = ctrnn;
    }

    private Braitenberg(sensor_1:number,sensor_2:number){
        //ASSUMPTION
        //assume sensors directly drive the opposite side's motor
        //a la braitenberg's robot
        var force_multiplier = 1;//0.6;
        var mu = 0.025; var sigma = 0.12;
        var noise = Rand(0,1) * sigma + mu
        var l_accelleration = (sensor_1 + noise) * force_multiplier;
        this.LeftMotorSpeed += l_accelleration;
        noise = Rand(0,1) * sigma + mu
        var r_accelleration = (sensor_2 + noise) * force_multiplier;
        this.RightMotorSpeed += r_accelleration;

        ///Friction application
        var friction = 0.1;
        this.LeftMotorSpeed = Clamp(0,10,this.LeftMotorSpeed-=friction);
        this.RightMotorSpeed = Clamp(0,10,this.RightMotorSpeed-=friction);
        //**** debug
        // this.RightMotorSpeed = 1.0;
        // this.LeftMotorSpeed = 1.1;
        //****
        // console.log('motor speeds ',this.LeftMotorSpeed,this.RightMotorSpeed);
    }

    private ComputeCtrnn(sensor_1:number,sensor_2:number){
        var mu = 0.25; var sigma = 0.25;
        var l_sensor = (sensor_1 + RandNoise(mu,sigma)) * this.LeftSensorGain;
        var r_sensor = (sensor_2 + RandNoise(mu,sigma)) * this.RightSensorGain;
        // console.log('l_sensor,r_sensor',l_sensor,r_sensor);

        //Set input
        this.Ctrnn.I[1] = l_sensor;
        this.Ctrnn.I[3] = r_sensor;

        var k = 5; //Evaluate CTRNN k times before taking the output
        var homeostatic_neuronal_ratio = 0;
        for(var i=0;i<k;i++){
            homeostatic_neuronal_ratio += this.Ctrnn.Evaluate();
        }

        this.HomeostaticNeuronalRatio += homeostatic_neuronal_ratio/k;
        this.CurrentPlasticNeuronalRatio = (1 - homeostatic_neuronal_ratio/k);
        //Get output and set to motors
        var q_left = LinScale(-1,1, this.Ctrnn.Z[5]); //remap to -1 to +1
        var q_right = LinScale(-1,1, this.Ctrnn.Z[7]); //remap to -1 to +1

        // console.log('q_left',q_left,'q_right',q_right);

        mu = 0; sigma = 0.25
        this.LeftMotorSpeed  = (q_left  + RandNoise(mu,sigma)) * this.LeftMotorGain;
        this.RightMotorSpeed = (q_right + RandNoise(mu,sigma)) * this.RightMotorGain;
        // console.log('Y',this.Ctrnn.Y)
    }

    private ComputeMotorSpeeds(sensor_1,sensor_2){
        // this.Braitenberg(sensor_1,sensor_2);
        this.ComputeCtrnn(sensor_1,sensor_2);
    }

    private ComputeIncidentAngle(
        base_point:Vector2D, //base point
        incident_point:Vector2D, //incident point
        light_source:Vector2D //light source
        ){
        var v1 = incident_point.sub(base_point);//normal to the tangent
        var v2 = light_source.sub(incident_point);//light ray incident vector
        var angle_rad = v1.angle(v2);
        return angle_rad;
    }

    private ComputeSensorResponse(v_solar:Vector2D,angle:number):number{
        //compute euclidean distance to the light source from sensors
        //square of which is inversely proportional to the sensor
        //activations

        //compute location on the sensors
        var a1 = angle + this.Theta;
        var xx1 = this.Position.x + this.Radius * Math.cos(ToRadians(a1));
        var yy1 = this.Position.y + this.Radius * Math.sin(ToRadians(a1));
        // var d1 = (xx1-x) * (xx1-x) + (yy1-y) * (yy1-y);
        // console.log('angle,dist',a1,d1)
        // console.log('xx yy',xx1,yy1);
        //---------

        var sensor = new Vector2D(xx1,yy1);

        // console.log('base',this.X,this.Y,'solar',x,y,'sensor',xx1,yy1)
        var incident_angle_rad = this.ComputeIncidentAngle(this.Position,sensor,v_solar);
        // console.log('incident_angle',ToDegrees(incident_angle_rad),incident_angle_rad);
        if (incident_angle_rad> ToRadians(-90) && incident_angle_rad < ToRadians(90)){
            var distance = sensor.sub(v_solar).norm();
            var d1 = distance * distance;
            //light energy on sensors obey the inverse square law
            var perceived_intensity1 = 2000000/d1;
            //console.log('intensity',perceived_intensity1)
            return perceived_intensity1;
        }
        else
            return 0;
        //---------
    }

    private ComputeSolarVicinityTime(v_solar:Vector2D,r_solar:number){
        // if this in this timestep, if the bot is closer to the source
        // more than 4* radius of the bot, it is considered as SolarVicinityTime
        var dist = v_solar.sub(this.Position);
        var dist_n = dist.norm();
        if(dist_n < r_solar * 4)
            this.SolarVicinityTime +=1.0;

        this.IntegralSolarDistanceOverTime += (1 / (1 + dist_n));
    }
    //update internal state
    //x & y are locations of the light source
    public Update(v_solar:Vector2D,r_solar:number){
        //compute euclidean distance to the light source from sensors
        //square of which is inversely proportional to the sensor
        //activations

        //compute location on the sensors
        var sensor_1 = this.ComputeSensorResponse(v_solar,this.SensorAngle);
        var sensor_2 = this.ComputeSensorResponse(v_solar,360-this.SensorAngle);
        //console.log('sensor output', sensor_1,sensor_2);

        //Feed the sensor inputs to the CTRNN and compute motor outputs
        if(this.IsInverted)
            this.ComputeMotorSpeeds(sensor_2,sensor_1);
        else
            this.ComputeMotorSpeeds(sensor_1,sensor_2);
        // console.log('motors left,right',this.LeftMotorSpeed,this.RightMotorSpeed);

        // var forward_speed = (this.LeftMotorSpeed+this.RightMotorSpeed) /2;
        // var rad_angle = (this.LeftMotorSpeed-this.RightMotorSpeed) / (2 * this.Radius);
        // var angle = ToDegrees(rad_angle);
        // console.log('angle, fwd speed',angle,forward_speed);
        // this.Theta += angle;

        // var a = ToRadians(this.Theta);
        // var xx = forward_speed * Math.cos(a);
        // var yy = forward_speed * Math.sin(a);
        // // console.log('xx,yy',xx,yy)
        // this.X += xx; this.Y += yy;

        //left motor displacement
        var rad_cw_delta_angle = this.LeftMotorSpeed / (2*this.Radius);
        var x_disp_l = this.LeftMotorSpeed / 2;

        var rad_ccw_detla_angle = this.RightMotorSpeed / (2*this.Radius);
        var x_disp_r = this.RightMotorSpeed /2;

        var total_angle_rad = rad_cw_delta_angle - rad_ccw_detla_angle;
        var total_angle = ToDegrees(total_angle_rad);

        //console.log('diff-angles', total_angle_rad, total_angle);
        var a = ToRadians(this.Theta);
        var xx = (x_disp_l + x_disp_r) * Math.cos(a);
        var yy = (x_disp_l + x_disp_r) * Math.sin(a);
        this.Theta += total_angle;

        //console.log('xx,yy',xx,yy)
        this.Position.x += xx; this.Position.y += yy;

        //Update metrics
        this.ComputeSolarVicinityTime(v_solar,r_solar);
        // console.log('update',this.Theta,this.X,this.Y);
    }

    public Draw(ctx:CanvasRenderingContext2D){
        var x = this.Position.x;
        var y = this.Position.y;

        ctx.save();
        ctx.translate(x,y);
        var angle_rad = ToRadians(this.Theta);
        ctx.rotate(angle_rad);
        ctx.translate(-x, -y);

        //draw circular body
        ctx.beginPath();
        ctx.arc(x,y,this.Radius,0,2.0*Math.PI);
        ctx.closePath();
        ctx.fillStyle = this.CurrentPlasticNeuronalRatio > 0.5 ? '#43C6DB': '#616D7E';
        ctx.fill();
        // ctx.save();

        //draw wheels
        ctx.fillStyle = 'brown';
        var wheel_width = 5; var wheel_length = 20;
        var xx = x - wheel_length / 2.0;
        var yy = y - this.Radius - wheel_width;
        ctx.fillRect(xx, yy, wheel_length, wheel_width);
        yy += this.Radius * 2.0 + wheel_width;
        ctx.fillRect(xx, yy, wheel_length, wheel_width);
        // ctx.save();

        //draw orientation line
        ctx.lineWidth = 1;
        ctx.strokeStyle = 'black';
        xx = x - this.Radius; yy = y;
        ctx.lineTo(xx,yy);
        ctx.stroke();
        // ctx.save();

        //draw sensors
        var sensor_r = 5;
        xx = x + this.Radius * Math.cos(ToRadians(this.SensorAngle));
        yy = y + this.Radius * Math.sin(ToRadians(this.SensorAngle))
        ctx.beginPath();
        ctx.arc(xx,yy,sensor_r,0,2.0*Math.PI);
        ctx.closePath();
        ctx.fillStyle = 'yellow';
        ctx.fill();

        var yy1 = y - this.Radius * Math.sin(ToRadians(this.SensorAngle))
        ctx.beginPath();
        ctx.arc(xx,yy1,sensor_r,0,2.0*Math.PI);
        ctx.closePath();
        ctx.fillStyle = 'yellow';
        ctx.fill();
        // ctx.save();

        //draw id
        ctx.font = "50px serif";
        var str = ""+(this.Id-1000);
        ctx.fillText(str, x-10, y+10);

        // //draw guides
        // ctx.moveTo(x,y);
        // ctx.lineTo(xx,yy)
        // //ctx.moveTo(xx,yy);
        // ctx.restore();
        // ctx.lineTo(500,500);
        // ctx.stroke();

        // ctx.restore();
        // ctx.restore();
        // ctx.restore();

        ctx.restore();
    }
}

//----------------------------------------------------------------------------
// Light Source class
//----------------------------------------------------------------------------
class Light implements ISituatedObject{
    constructor(public Id:number, public Position:Vector2D, public Radius:number){
    }

    //update internal state
    public Update(vect:Vector2D){
    }

    public Draw(ctx:CanvasRenderingContext2D){
        ctx.beginPath();
        ctx.arc(this.Position.x,this.Position.y,this.Radius,0,2.0*Math.PI);

        // // flicker the light source intensity. small animation
        // // to show it is alive
        // var delta = parseInt(Rand(-8,8)+'');
        // var r = 255 + delta; var g = 69 + delta; var b = 0 + delta;
        // var color = RGB2Hex(r,g,b);
        // //console.log(color);
        var color = 'red';

        ctx.fillStyle = color;
        ctx.fill();
        ctx.lineWidth = 2;
        ctx.strokeStyle = 'black';
        ctx.stroke();
        ctx.closePath();

        //draw the vicinity
        ctx.beginPath();
        ctx.setLineDash([5]);
        ctx.arc(this.Position.x,this.Position.y,this.Radius*4,0,2.0*Math.PI);
        ctx.lineWidth = 2;
        ctx.strokeStyle = 'red';
        ctx.stroke();
        ctx.closePath();

    }
}

//----------------------------------------------------------------------------
// Genetic State Machine states
//----------------------------------------------------------------------------
enum EvolutionaryLifeCycle{
    Genesis, //bots get their chromosome set and get created
    EvaluateFitness, //bots live out in their designated environment
    Evolve, //evolve the next generation
    TestLongTermStability, //After K generations of evolutions, take the best bots and evaluate
                          //their performance
    TestInversionStability //After K generations of normal phototaxis, sensors are inverted and tested for
                            //phototaxis
}

//----------------------------------------------------------------------------
// Simulation Class
//----------------------------------------------------------------------------
class RobotSimulation{
    private CachedCanvas:CanvasRenderingContext2D;
    private SituatedObjects: Array<ISituatedObject> = new Array<ISituatedObject>();
    public FrameCounter:number = 0;
    private SimulationStepCounter = 0;
    private MaxSimulationSteps:number = 500; //Max Simulation steps per generation

    private LifeCycle:EvolutionaryLifeCycle = EvolutionaryLifeCycle.Genesis;
    private GA = new SteadyStateGeneticAlgorithm(60,0.7,0.7,148,64);
    private MSqFitnessChart:Chart=null;
    private SolarFlockingChart:Chart=null;
    constructor(public Width: number, public Height: number) {
        var s = this.LoadChromosomes();
        if(s.trim()!=""){
            console.log('loading Chromosomes.. ', s.length);
            this.GA.Deserialize(s);
            var solar_pos = new Vector2D(Width/2,Height/2);
            var solar = new Light(1000,solar_pos,50);
            this.CreateBots(solar);
            this.LifeCycle = EvolutionaryLifeCycle.TestLongTermStability;
            console.log('loaded');
        }
    }

    private MapChromosomeToNetwork(chromo:Array<number>){
        //ChromosomeMap
        //Weights(64),PlasticEta(64),Biases(8),Tau(8),LMotorGain(1),RMotorGain(1),LSensorGain(1),RSensorGain(1)

        var n = 8;// number of neurons
        var kk = 0;

        //weight matrix
        var w = new Array<Array<number>>(n);
        for(var i=0;i<n;i++){
            w[i] = new Array<number>(n);
            for(var j=0;j<n;j++){
                w[i][j] = LinScale(-8,8,chromo[kk++]); //remap range[0,1] to [-8,8]
            }
        }
        //plastic eta matrix
        var eta = new Array<Array<number>>(n);
        for(var i=0;i<n;i++){
            eta[i] = new Array<number>(n);
            for(var j=0;j<n;j++){
                eta[i][j] = LinScale(-0.9,0.9,chromo[kk++]); //remap range[0,1] to [-0.9,0.9]
            }
        }

        var b = new Array<number>(n);
        var y = new Array<number>(n);
        var tau = new Array<number>(n);
        var input = new Array<number>(n);
        for(var i=0;i<n;i++){
            b[i] = LinScale(-3,3,chromo[kk+i]);//remap range[0,1] to [-3,3];
            tau[i] = LinScale(0.4,4,chromo[kk+n+i]); //remap range [0,1] to (0.4,4);
        }

        //Plastic Rules
        kk = 148; //offset to the plastic rule segment
        var rules= new Array<Array<number>>(n);
        for(var i=0;i<n;i++){
            rules[i] = new Array<number>(n);
            for(var j=0;j<n;j++){
                rules[i][j] = chromo[kk++];
            }
        }

        return  new Network(w,b,tau,rules,eta)//Create the ctrnn
    }

    public InitRender(ctx:CanvasRenderingContext2D){
        //create grid
        var d = 10;
        for(var x = 0; x <=this.Width; x += d){
				ctx.moveTo(x,0);
				ctx.lineTo(x,this.Height);
			}
			for(var y = 0; y <=this.Height; y+=d){
				ctx.moveTo(0,y);
				ctx.lineTo(this.Width,y);
			}
			ctx.strokeStyle = "gainsboro";
			ctx.lineWidth = 0.5;
			ctx.stroke();

        this.CachedCanvas = ctx;
    }

    private CreateBots(solar:Light){
        //read chromosome from GA and create bots
        this.SituatedObjects = new Array<ISituatedObject>();

        //create a light source
        //------------------------
        // var id= 1000;
        // // var x = Rand(0,this.Width);
        // // var y = Rand(0,this.Height);
        // var x_s = this.Width/2; var y_s = this.Height/2; var r = 20;
        // var solar = new Light(id++, x_s, y_s, r);
        this.SituatedObjects.push(solar);
        var id = solar.Id;
        //----------------------
        //ChromosomeMap
        //Weights(64),PlasticEta(64),Biases(8),Tau(8),LMotorGain(1),RMotorGain(1),LSensorGain(1),RSensorGain(1)
        var yy = 0;
        for(var chromo of this.GA.Population){
            // console.log('chromo',chromo)
            var r = 50;

            //random position whole space
            //---------------------------
            // var x = Rand(this.Width/4,3*this.Width/4);
            // var y = Rand(this.Height/4,3*this.Height/4);
            // var angle = Rand(-180,180); //facing angle

            //lined up on the left wall
            //-------------------------
            // var x = r+10;
            // var y = 10 + yy * r*2;
            // yy++;
            //var angle = 0;

            //around the solar
            //----------------
            var pos = RandomCartesianCoordinate(1000);
            var angle = Rand(0,360);
            pos = pos.add(new Vector2D(this.Width/2, this.Height/2));

            var ctrnn = this.MapChromosomeToNetwork(chromo);
            var robot = new Bot(id++,pos,r,angle,ctrnn);
            var k = 64+64+8+8;//add offset of weights(64),plasticEta(64), biases(8) and tau(8)
            robot.LeftMotorGain = ExpScale(0.1,10,chromo[k++]);     //remap to [0.1,10] exponentially scaled
            robot.RightMotorGain = ExpScale(0.1,10,chromo[k++]);    //remap to [0.1,10] exponentially scaled
            robot.LeftSensorGain = ExpScale(0.1,10,chromo[k++]);    //remap to [0.1,10] exponentially scaled
            robot.RightSensorGain = ExpScale(0.1,10,chromo[k++]);   //remap to [0.1,10] exponentially scaled

            robot.Di = pos.sub(solar.Position).norm(); //dist to solar
            this.SituatedObjects.push(robot);
        }
        //transit state
        if(this.GA.CurrentGeneration >600){
            this.LifeCycle = EvolutionaryLifeCycle.TestLongTermStability;
            this.SolarFlockingChart.Clear();
            //Save GA bits
            this.SaveChromosomes(this.GA.SerializeTo());
        }
        else
            this.LifeCycle = EvolutionaryLifeCycle.EvaluateFitness;
    }
    private SaveChromosomes(data:string){
        var list = document.getElementById('output');
        var entry = document.createElement('li');
        var input = document.createElement('textarea');
        input.readOnly = true;
        input.value = data;
        entry.appendChild(input);
        list.appendChild(entry);
    }
    private LoadChromosomes():string{
        var input = <HTMLTextAreaElement>document.getElementById('input-text');
        return input.value;
    }

    private GetSolarPosition():Vector2D{
        /*/
        var total = 0;
        for(var i=0;i<=this.CurrentSolarSegment;i++){
            total += this.SolarDurations[i];
        }
        if(this.SimulationStepCounter>=total){
            //compute fitness
            var v_solar = this.SolarPositions[this.CurrentSolarSegment]
            for(var i=1;i<this.SituatedObjects.length;i++){
                    var bot = <Bot>(this.SituatedObjects[i]); //robot
                    bot.Df = ((new Vector2D(bot.X,bot.Y)).sub(v_solar)).norm(); //dist to solar
            }

            this.CurrentSolarSegment++;
        }


        return this.SolarPositions[this.CurrentSolarSegment];
        /*/
        var solar = this.SituatedObjects[0]
        return solar.Position
        //*/
    }

    private EvaluateBots(ctx:CanvasRenderingContext2D, canvas:HTMLCanvasElement){
        ctx.save();
        ctx.clearRect(0,0,this.Width,this.Height);
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.drawImage(this.CachedCanvas.canvas, 0, 0);

        var solar = <Light>this.SituatedObjects[0]
        for (var so of this.SituatedObjects) {
                so.Update(solar.Position,solar.Radius);//Update the bots
                so.Draw(ctx);//Draw the bots
        }
        ctx.restore();

        //Show the bread crumbs of the bots
        var ctx_ = this.CachedCanvas;
        for(var i=1;i<this.SituatedObjects.length;i++){
                var robot =this.SituatedObjects[i] //robot
                ctx_.fillStyle = 'lightgrey';
                ctx_.fillRect(robot.Position.x, robot.Position.y, 3, 3);
        }

        //Compute metrics and transit state
        //---------------------------------
        this.SimulationStepCounter++;
        if(this.SimulationStepCounter > this.MaxSimulationSteps){//simulation finished
            this.SimulationStepCounter=0; //reset counter
            this.LifeCycle = EvolutionaryLifeCycle.Evolve;//transit state

            //record final dist
            for(var i=1;i<this.SituatedObjects.length;i++){
                    var bot = <Bot>(this.SituatedObjects[i]); //robot
                    bot.Df = bot.Position.sub(solar.Position).norm(); //dist to solar
            }
        }
    }

    private LTStabilityTestGenerationCount:number = 0;
    private TestLTStability(ctx:CanvasRenderingContext2D, canvas:HTMLCanvasElement){
        //this would be called on every frame
        ctx.save();
        ctx.clearRect(0,0,this.Width,this.Height);
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.drawImage(this.CachedCanvas.canvas, 0, 0);

        var solar = <Light>this.SituatedObjects[0]
        for (var so of this.SituatedObjects) {
                so.Update(solar.Position,solar.Radius);//Update the bots
                so.Draw(ctx);//Draw the bots
        }
        ctx.restore();

        this.SimulationStepCounter++;
        if(this.SimulationStepCounter > this.MaxSimulationSteps){//simulation finished
            this.SimulationStepCounter=0; //reset counter
            this.LTStabilityTestGenerationCount++;

            this.ComputeVicinityMetric();

            //Change Solar location
            var replaced =false;
            while(!replaced){
                var disp = RandomCartesianCoordinate(RandomInt(0,500));
                var new_pos = solar.Position.add(disp);
                if(new_pos.x > 0 && new_pos.x < this.Width
                            && new_pos.y > 0  && new_pos.y < this.Height){
                    solar.Position = new_pos;
                    replaced = true;
                }
            }
        }

        if(this.LTStabilityTestGenerationCount>100){
            this.LTStabilityTestGenerationCount = 0;

            //transition to the next stage
            //States transition back and forth between normal and inverted sensors
            if(this.LifeCycle == EvolutionaryLifeCycle.TestInversionStability){
                this.LifeCycle = EvolutionaryLifeCycle.TestLongTermStability;
                console.log('Sensors Normal')
                for(var b of this.SituatedObjects.slice(1)){
                    var bot = <Bot>b;
                    bot.IsInverted = false;
                }
            }else{
                this.LifeCycle = EvolutionaryLifeCycle.TestInversionStability;
                console.log('Sensors Inverted')
                for(var b of this.SituatedObjects.slice(1)){
                    var bot = <Bot>b;
                    bot.IsInverted = true;
                }
            }
        }
    }

    private ComputeFitnessMetrics(){
        var bot_count = this.SituatedObjects.length-1;
        var fitness = new Array<number>(bot_count); //create fitness vector
        var fit_total = 0;

        var solar = <Light>(this.SituatedObjects[0]);
        var solar_flock_percent = 0;
        for(var i=1;i<bot_count+1;i++){
            var bot = <Bot>(this.SituatedObjects[i])

            var Fd = 1 - bot.Df/bot.Di
            // fitness[i-1] = (bot.SolarVicinityTime + bot.IntegralSolarDistanceOverTime);
            var Fp = bot.SolarVicinityTime / this.MaxSimulationSteps;
            var Fh = bot.HomeostaticNeuronalRatio / this.MaxSimulationSteps;
            bot.HomeostaticNeuronalRatio =0; //reset

            // var W1 = 0.21; //Weights add up to 1
            // var W2 = 0.64; //0.64 - 0.68 range
            // var W3 = 0.15; //0.15 - 0.2 range
            var W1 = 0.2; //Weights add up to 1
            var W2 = 0.3; //0.64 - 0.68 range
            var W3 = 0.5; //0.15 - 0.2 range

            // var f = (W1*Fd + W2*Fp + W3*Fh);
            var f = W1*Fd+ 1000 * (W2*Fp * W3*Fh);
            fitness[i-1] = f;
            // console.log(i,'fitness,Fd',f);

            fit_total += (f*f);

            if(solar.Position.sub(bot.Position).norm()< solar.Radius*4){
                solar_flock_percent++;
                console.log('flocked bot',bot.Id)
            }
        }

        //metrics display
        var avg_sq_fitness = fit_total/bot_count;
        this.MSqFitnessChart.Update(avg_sq_fitness);
        console.log('mean sq fitness',this.GA.CurrentGeneration,',',avg_sq_fitness);

        solar_flock_percent = solar_flock_percent/bot_count;
        this.SolarFlockingChart.Update(solar_flock_percent);
        console.log('solar flock percent',this.GA.CurrentGeneration,',',solar_flock_percent);
        return fitness;
    }

    private ComputeVicinityMetric(){
        var bot_count = this.SituatedObjects.length-1;

        var solar = <Light>(this.SituatedObjects[0]);
        var solar_flock_percent = 0;
        var homeostatic_percent = 0.0;
        var bot_list="";
        for(var i=1;i<bot_count+1;i++){
            var bot = <Bot>(this.SituatedObjects[i])
            var dist = solar.Position.sub(bot.Position).norm()
            if(dist < solar.Radius*4){
                solar_flock_percent++;
                // console.log('flocked bot',bot.Id)
                bot_list = bot_list + " "+ bot.Id;
            }

            var Fh = bot.HomeostaticNeuronalRatio / this.MaxSimulationSteps;
            console.log(bot.Id,'Fh', Fh,'dist', dist);
            bot.HomeostaticNeuronalRatio = 0; //reset -NOTE this is bad, resetting here
                                              //(only because this is not used anywhere else)
            homeostatic_percent += Fh; //(approxEqual(Fh,0) ? 0: 1);
        }

        var avg_homeostatic_bots = homeostatic_percent / bot_count;
        //metrics display
        solar_flock_percent = solar_flock_percent/bot_count;
        this.SolarFlockingChart.Update(solar_flock_percent);
        console.log(this.LTStabilityTestGenerationCount,',',
            'solar-flock-percent',
            solar_flock_percent,
            avg_homeostatic_bots,
            bot_list
            );
    }

    private Evolve(){
        var fitness = this.ComputeFitnessMetrics();
        this.GA.Evolve(fitness);
        this.LifeCycle = EvolutionaryLifeCycle.Genesis;//transit state
    }

    private SolarDurations:Array<number>=null;
    private SolarPositions:Array<Vector2D>=null;
    private CurrentSolarSegment = 0;
    private MaxSolarSegments = 5;
    private CreateSolarSchedule(){
        //initialize if not already done
        if(this.SolarDurations==null){
            this.SolarPositions = new Array<Vector2D>(this.MaxSolarSegments);
            this.SolarDurations = new Array<number>(this.MaxSolarSegments);
        }

        var x_s = this.Width/2; var y_s = this.Height/2;
        for(var i=0;i<this.MaxSolarSegments;i++){
            var T = this.MaxSimulationSteps/this.MaxSolarSegments;
            this.SolarDurations[i] = RandomInt(0.75*T,1.25*T);
            x_s += RandomInt(0,2000);
            y_s += RandomInt(-2000,2000);
            this.SolarPositions[i] = new Vector2D(x_s,y_s);
        }
    }

    private SolarDelayUponInversion_Counter = 0;

    public Update(ctx:CanvasRenderingContext2D, canvas:HTMLCanvasElement){
        var solar_pos = new Vector2D(this.Width/2,this.Height/2)
        //experiment:
        //----------------------------------------------------------------
        //After K generations we move the solar randomly within some range
        if(this.GA.CurrentGeneration>50){
            var new_pos = RandomCartesianCoordinate(RandomInt(0,500));
            solar_pos = solar_pos.add(new_pos);
        }
        //----------------------------------------------------------------
        var solar = new Light(1000,solar_pos,50);

        switch (this.LifeCycle) {
            case EvolutionaryLifeCycle.Genesis:{//single timestep
                this.CreateBots(solar);
                // this.CreateSolarSchedule();
            }break;
            case EvolutionaryLifeCycle.EvaluateFitness:{//multiple timesteps
                this.EvaluateBots(ctx,canvas);
            }break;
            case EvolutionaryLifeCycle.Evolve:{//single timestep
                this.Evolve();
            }break;
            case EvolutionaryLifeCycle.TestInversionStability:
            case EvolutionaryLifeCycle.TestLongTermStability:{
                this.TestLTStability(ctx,canvas);
            }break;
            // case EvolutionaryLifeCycle.TestInversionStability:{//multiple timesteps{//multiple timesteps
            //     //add a delay before solar appear  -FIX THIS, bug - bots are not simulated during
            //     //                                  this time.. which they should be and
            //     //                                  solar should be at infinity
            //     if(this.SolarDelayUponInversion_Counter> 200){
            //         this.SolarDelayUponInversion_Counter = 0;
            //         this.LifeCycle = EvolutionaryLifeCycle.TestLongTermStability;
            //     }else{
            //         this.SolarDelayUponInversion_Counter++;
            //     }
            // }break;
            default:
                break;
        }
        //draw metrics
        ctx.font = "48px serif";
        var str = "SimulationStepCounter = " + this.SimulationStepCounter
                    + " Gen Count = " + this.GA.CurrentGeneration
                    + " Max Fitness = " + (Math.max.apply(Math, this.GA.Fitness)).toFixed(3);
        ctx.fillText(str, 10, 50);

        //create chart
        if(this.MSqFitnessChart==null)
            this.MSqFitnessChart = new Chart(ctx,new Vector2D(2500,10),new Vector2D(640,240),ChartType.Log);
        else
            this.MSqFitnessChart.Draw();

        if(this.SolarFlockingChart==null)
            this.SolarFlockingChart = new Chart(ctx,new Vector2D(2500+640+100,10), new Vector2D(640,240),ChartType.Linear);
        else
            this.SolarFlockingChart.Draw();

        this.FrameCounter += 1;
    }
}

//----------------------------------------------------------------------------
// Main hook
//----------------------------------------------------------------------------
function exec() {
    var canv = document.createElement("canvas");

    var canvasTemp = document.createElement("canvas");
    var CanvasWidth = 4000; var CanvasHeight = 4000;
    canv.width = canvasTemp.width = CanvasWidth;
    canv.height = canvasTemp.height = CanvasHeight;
    document.body.appendChild(canv);
    var ctx = canv.getContext("2d");
    var tctx = canvasTemp.getContext("2d");
    var sim = new RobotSimulation(CanvasWidth, CanvasHeight);
    sim.InitRender(tctx); // cached rendering done here
    //return antsim.render(defaultScene(), ctx, CanvasWidth, CanvasHeight);
    setInterval(() => sim.Update(ctx,canv), 20);
}
//----------------------------------------------------------------------------
