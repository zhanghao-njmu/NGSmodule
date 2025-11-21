/**
 * Task List Page - Task monitoring with real-time updates
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Button,
  Table,
  Tag,
  Progress,
  Space,
  Select,
  Statistic,
  Row,
  Col,
} from 'antd'
import {
  PlusOutlined,
  SyncOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ClockCircleOutlined,
  PlayCircleOutlined,
  StopOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useTaskStore } from '../../store/taskStore'
import { useProjectStore } from '../../store/projectStore'
import { websocketService } from '../../services/websocket.service'
import type { Task, TaskStatus } from '../../types/task'
import dayjs from 'dayjs'

const { Option } = Select

export const TaskList: React.FC = () => {
  const {
    tasks,
    stats,
    loading,
    fetchTasks,
    fetchStats,
    cancelTask,
    handleWebSocketMessage,
  } = useTaskStore()
  const { projects, fetchProjects } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')

  useEffect(() => {
    fetchProjects()
    fetchStats()

    // Setup WebSocket for real-time updates
    const token = localStorage.getItem('auth_token')
    if (token) {
      websocketService.connect(token)
      websocketService.addMessageHandler(handleWebSocketMessage)
    }

    return () => {
      websocketService.disconnect()
    }
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchTasks({ project_id: selectedProject })
    } else {
      fetchTasks()
    }
  }, [selectedProject])

  const getStatusConfig = (status: TaskStatus) => {
    const configs = {
      pending: { color: 'default', icon: <ClockCircleOutlined /> },
      running: { color: 'processing', icon: <SyncOutlined spin /> },
      completed: { color: 'success', icon: <CheckCircleOutlined /> },
      failed: { color: 'error', icon: <CloseCircleOutlined /> },
      cancelled: { color: 'warning', icon: <StopOutlined /> },
    }
    return configs[status] || configs.pending
  }

  const columns: ColumnsType<Task> = [
    {
      title: 'Task Name',
      dataIndex: 'task_name',
      key: 'task_name',
      render: (name) => <span style={{ fontWeight: 500 }}>{name}</span>,
    },
    {
      title: 'Type',
      dataIndex: 'task_type',
      key: 'task_type',
      width: 120,
      render: (type) => type && <Tag color="blue">{type}</Tag>,
    },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      width: 120,
      render: (status: TaskStatus) => {
        const config = getStatusConfig(status)
        return (
          <Tag color={config.color} icon={config.icon}>
            {status.toUpperCase()}
          </Tag>
        )
      },
    },
    {
      title: 'Progress',
      dataIndex: 'progress',
      key: 'progress',
      width: 150,
      render: (progress) => (
        <Progress
          percent={Math.round(progress)}
          size="small"
          status={progress === 100 ? 'success' : 'active'}
        />
      ),
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 100,
      render: (_, record) => (
        <Space>
          {record.status === 'running' && (
            <Button
              size="small"
              danger
              icon={<StopOutlined />}
              onClick={() => cancelTask(record.id)}
            >
              Cancel
            </Button>
          )}
        </Space>
      ),
    },
  ]

  return (
    <div>
      {/* Statistics */}
      <Row gutter={16} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={8} lg={4}>
          <Card>
            <Statistic
              title="Total"
              value={stats?.total_tasks || 0}
              valueStyle={{ color: '#1890ff' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={8} lg={4}>
          <Card>
            <Statistic
              title="Running"
              value={stats?.running_tasks || 0}
              valueStyle={{ color: '#fa8c16' }}
              prefix={<SyncOutlined spin />}
            />
          </Card>
        </Col>
        <Col xs={24} sm={8} lg={4}>
          <Card>
            <Statistic
              title="Completed"
              value={stats?.completed_tasks || 0}
              valueStyle={{ color: '#52c41a' }}
              prefix={<CheckCircleOutlined />}
            />
          </Card>
        </Col>
        <Col xs={24} sm={8} lg={4}>
          <Card>
            <Statistic
              title="Failed"
              value={stats?.failed_tasks || 0}
              valueStyle={{ color: '#ff4d4f' }}
              prefix={<CloseCircleOutlined />}
            />
          </Card>
        </Col>
      </Row>

      <Card style={{ marginBottom: 16 }}>
        <Space style={{ width: '100%', justifyContent: 'space-between' }}>
          <Select
            placeholder="All Projects"
            style={{ width: 300 }}
            value={selectedProject || undefined}
            onChange={setSelectedProject}
            allowClear
          >
            {projects.map((p) => (
              <Option key={p.id} value={p.id}>
                {p.name}
              </Option>
            ))}
          </Select>
          <Button type="primary" icon={<PlusOutlined />}>
            New Task
          </Button>
        </Space>
      </Card>

      <Card>
        <Table
          columns={columns}
          dataSource={tasks}
          rowKey="id"
          loading={loading}
          pagination={{ pageSize: 20 }}
        />
      </Card>
    </div>
  )
}

export default TaskList
